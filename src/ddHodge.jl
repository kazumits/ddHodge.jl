module ddHodge

export schurselect, basischange, kNNGraph, ddHodgeWorkflow

using LinearAlgebra, SparseArrays, Graphs
using KrylovKit: linsolve, eigsolve
using NearestNeighbors: KDTree, knn

module DEC
    function delaunayGraph end
    function stardiags end
end

module GPU
    # Placeholder of `CUDA.CUSPARSE.CuSparseMatrixCSR`
    function cuspm end
    # Placeholder of `CUDA.CuVector`
    function cuvec end
end

"""
    graphgrad(g)

Construct grad operator from graph `g`

The direction is implicitly determined by node order, i.e., `i -> j` if `i < j`.

# Example
```jldoctest
julia> using Graphs

julia> g = complete_graph(3)
{3, 3} undirected simple Int64 graph

julia> collect(edges(g))
3-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
 Edge 1 => 3
 Edge 2 => 3

julia> d0 = ddHodge.graphgrad(g)
3×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 6 stored entries:
 -1.0   1.0   ⋅ 
 -1.0    ⋅   1.0
   ⋅   -1.0  1.0
```
"""
function graphgrad(g::SimpleGraph{Int})
    I = sparse(1:ne(g), src.(edges(g)), -1.0, ne(g), nv(g))
    J = sparse(1:ne(g), dst.(edges(g)),  1.0, ne(g), nv(g))
    I + J
end

# Approximate the integration of flux along edges
function calcwt(g::SimpleGraph{Int}, pts::AbstractMatrix, vel::AbstractMatrix)
    map(edges(g)) do e
        i, j = src(e), dst(e)
        dx = pts[:,j] - pts[:,i] # i to j
        (vel[:,i] + vel[:,j])'*dx/2.0
    end
end

"""
    kNNGraph(X::Matrix, k::Int) -> (g::SimpleGraph, kdtree::KDTree)

Graph construction using k-NN
"""
function kNNGraph(X::AbstractMatrix, k::Int)
    N = size(X,2)
    kdtree = KDTree(X)
    idxs, dists = knn(kdtree, X, k+1, true)
    g = SimpleGraph(N)
    for (i,v) in enumerate(idxs)
        for j in 2:(k+1)
            # disconnect if zero distance
            if dists[i][j] > 0.0
                add_edge!(g,v[begin],v[j])
            end
        end
    end
    #rem_vertices!(g,findall(degree(g).==0))
    #@assert is_connected(g) # stop when g is not connected
    g, kdtree, idxs, dists
end
 
# Coefficients of directional derivatives w.r.t. Hessian
function chesse(v::Vector{T}) where T<:AbstractFloat
    n = length(v)
    c = fill(0.0::T, Int(n*(n-1)/2))
    k = 1
    for i in 1:n-1, j in i+1:n
        c[k] = 2*v[i]*v[j]
        k+=1
    end
    [v.^2; c]
end

# Preparing the coefficients to solve batch Hessian estimation
function localHessNS(
    g::SimpleGraph{Int},
    i::Int,
    pts::AbstractMatrix,
    vel::AbstractMatrix;
    B=I::Matrix{Float64}
)
    nei = neighbors(g,i)
    if length(nei) == 0
        return( # values for orphan node
            E = zeros(size(B,2)),
            y = 0.0, d = 0.0
        )
    end
    ηs = [pts[:,j] - pts[:,i] for j in nei]
    d = norm.(ηs)
    #E = B'*mapfoldl(normalize,hcat,ηs)
    E = B'*hcat(normalize.(ηs)...)
    # vi - vj since -gradU is the flow
    y = [η'*(vel[:,i]-vel[:,j]) for (η, j) in zip(ηs, nei)]./d.^2
    (E=E, y=y, d=d)
end

# Reshape the result into the Hessian matrix form
function chesse2mat(x::AbstractVector, d::Int)
    H = fill(0.0, d, d)
    H[diagind(H)] = x[1:d]
    k = d+1
    for i in 1:d-1, j in i+1:d
        H[i,j] = x[k]
        #H[j,i] = x[k] # no need for Symm.
        k += 1
    end
    Symmetric(H)
end

# Get dual velocities in the given span
function dualvels(V::AbstractMatrix, spans::Vector{<:AbstractMatrix}; d=[1,2]::Vector{Int})
    hcat(map(1:length(spans)) do i
        S = spans[i][:,d]
        S*[0 -1; 1 0]*S'*V[:,i]
    end...)
end

# General purpose function:
#  Propagate values to neighboring vertices from root (along BFS tree)
function propagate(f::Function, g::AbstractGraph, init::T; root::Int=1) where T
    #st = Array{typeof(init),1}(undef,nv(g))
    st = Array{T}(undef,nv(g))
    st[root] = init
    tree = bfs_tree(SimpleGraph(g),root) # bfs is better than dfs?
    parents = [root]
    for p in parents
        children = neighbors(tree,p)
        for c in children
            st[c] = f(p,c,st[p]) # assign some f(p,c,v)
            push!(parents,c)
        end
    end
    return st
end

# Get subspace of a i-th point by local PCA
function getspan(g::SimpleGraph{Int}, pts::AbstractMatrix, i::Int, d::Int; ref=nothing)
    nei = neighbors(g,i)
    nnei = length(nei)
    if nnei == 0
        return zeros(size(pts,1),1)
    end
    ηs = [pts[:,j] - pts[:,i] for j in nei]
    # hack for zero distance
    dists = norm.(ηs)
    #tooclose = dists .< 1e-16
    #nei = nei[.!tooclose]
    #ηs = ηs[.!tooclose]
    #dists = dists[.!tooclose]
    # higher weight for short distance
    wts = dists.^-2
    U = svd(Matrix(hcat(wts.*ηs...))).U
    # concatenate zeros if deficient
    U = size(U,2) < d ? hcat(U,zeros(size(pts,1),d-size(U,2) )) : U[:,1:d]
    if !isnothing(ref)
        # check whether Tᵢⱼ: Tᵢ ⟶ Tⱼ in SO(n)
        # note that sign(det(UΣVᵗ)) = sign(det(UVᵗ)) since detΣ>0
        if det(U'*ref) < 0
            # flipping to force the consistent orientations.
            U[:,end] = -U[:,end]
        end
    end
    U
end

# Get all spans while aligning orientations of tangent spaces
function oriented_spans(g::SimpleGraph{Int}, pts::AbstractMatrix, d::Int, root::Int)
    fdim, N = size(pts)
    rspan = getspan(g,pts,root,d,ref=Matrix(I,fdim,d))
    spans = propagate((pa,ch,vpa) -> getspan(g,pts,ch,d,ref=vpa), g, rspan, root=root)
    for i in 1:N
        if !isassigned(spans,i) # span of zero-degree vertex
            spans[i] = Matrix(1.0I, fdim, d) # any full rank matrix
        end
    end
    spans
end

# Storing the restriction map
struct Rmap 
    Oi ::Matrix{Float64} # transfer i -> e 
    Oj ::Matrix{Float64} # transfer j -> e
    #Ui ::Matrix{Float64} # span at i
    #Uj ::Matrix{Float64} # span at j
    pa ::Array{Float64}  # principal angles
end

# Calculate subspace alignment
# [ref] Singer, A. & Wu, H.-T. Vector diffusion maps and the connection Laplacian. Communications on Pure and Applied Mathematics 65, 1067–1144 (2012).
function connectspans(Ui::AbstractMatrix, Uj::AbstractMatrix)
    sv = svd(Matrix(Ui'Uj))
    pa = acos.(clamp.(sv.S,0,1)) # Principal angles
    Rmap(sv.U', sv.Vt, pa) # Oi=U', Oj=V' since Oi'*Oj=U*V'
end

# Gradient operater of sheaf over graph
# Cellular sheaf over graph, Hansen & Ghrist (2020)
# todo: applicability for variable sheaf dimensions
# remark: valid under the assumption of src(e) < dst(e)
function d0sheaf(g::SimpleGraph{Int}, rmaps::Dict{Tuple{Int,Int},Rmap})
    d  = length(first(rmaps)[2].pa)
    Is = Vector{Int}(undef,0)
    Js = Vector{Int}(undef,0)
    Vs = Vector{Float64}(undef,0)
    sizehint!.([Is,Js,Vs], 2*d*d*ne(g))
    for (i,e) in enumerate(edges(g))
        Fue = rmaps[Tuple(e)].Oi
        Fve = rmaps[Tuple(e)].Oj
        append!(Is, repeat((1:d) .+ d*(i-1), 2d))
        append!(
            Js,
            repeat((1:d) .+ d*(src(e)-1), inner=d), # Fue
            repeat((1:d) .+ d*(dst(e)-1), inner=d)  # Fve
        )
        append!(Vs, vec(Fue), vec(-Fve))
    end
    sparse(Is, Js, Vs, d*ne(g), d*nv(g))
end

"""
    alignSpaces(graph, subspaces) -> (alignedSpaces, sheafGrad, restrictionMaps)

Solves subspace alignment problem by minimizing Rayleigh quotient of sheaf (connection) laplacian ``\\Delta_s``

```math
\\min_{\\bm x} \\sum_{(u,v) \\in E} \\|O_u \\bm x_u - O_v \\bm x_v\\|^2 = \\bm x^\\top \\Delta_s \\bm x
```
where ``O_u, O_v`` are the restriction maps of parallel transport matrix.

"""
function alignSpaces(
    g::SimpleGraph{Int}, ss::Vector{<:AbstractMatrix};
    krytol=1e-32::Float64, useCUDA=false::Bool
)
    N = nv(g)
    damb, d = size(ss[begin])
    rms = Dict([Tuple(e) => connectspans(ss[src(e)],ss[dst(e)]) for e in edges(g)])
    d0 = d0sheaf(g,rms) # Sheaf gradient
    cuspm = useCUDA ? GPU.cuspm : identity
    cuvec = useCUDA ? GPU.cuvec : identity
    L = cuspm(d0'd0)
    #  x0 = cuvec(rand(N*d))
    x0 = cuvec(fill(1.0,N*d))
    #if usenormalL # faster calculation for Hessian?
        #Dhi = spdiagm(repeat(degree(g).^(-1/2),inner=d))
        #L = Dhi*L*Dhi
    #end
    vals, vecs, = eigsolve(L,x0,d,:SR,issymmetric=true,tol=krytol)
    vecs = collect.(vecs) # to cpu
    vec2dir(x) = normalize.(eachcol(reshape(x,(d,N))))
    # Parallel tangent vectors
    tvs = [hcat(ss.*vec2dir(v)...) for v in vecs[1:d]]
    # Aligned subspaces (pick each tangent vecs)
    #  sa = [hcat([tvs[i][:,j] for i in 1:d]...) for j in 1:N]
    sa = collect(eachslice(reshape(vcat(tvs...),damb,d,N),dims=3))
    sa, d0, rms
end

# Batch Hessian matrix estimation
function batchHess(
    g::SimpleGraph{Int}, pts::AbstractMatrix, vel::AbstractMatrix, Δ₀::SparseMatrixCSC;
    d=size(pts,1)::Int, λ=0.1::AbstractFloat, Ψ=fill(I,nv(g))::Array{AbstractMatrix},
    maxiter=5000::Int, useCUDA=false::Bool
)
    N = nv(g)
    # Prepare matrices for batch solve
    lh = [localHessNS(g,i,pts,vel,B=Ψ[i]) for i in 1:N]
    # D = I
    #y = mapfoldl(h -> h.y, vcat, lh)
    y = vcat(get.(lh,:y,missing)...)
    Ms = [sparse(mapslices(chesse,h.E,dims=1)') for h in lh]
    X = blockdiag(Ms...)
    #X = mapfoldl(h -> sparse(mapslices(chesse,h.E,dims=1)'), blockdiag, lh)
    p = size(Ms[begin],2)
    # L = -Δ₀ since (∇β)'(∇β) = -β'Δ₀β 
    L = kron(-Δ₀,sparse(I,p,p)) # shape of Hessians should be same
    #sol, = linsolve(X'D*X + λ*L, X'D*y, isposdef=true)
    cuspm = useCUDA ? GPU.cuspm : identity
    cuvec = useCUDA ? GPU.cuvec : identity
    A = cuspm(X'X + λ*L)
    b = cuvec(X'y)
    x0 = cuvec(rand(length(b)))
    sol, stats = linsolve(A,b,x0,isposdef=true,maxiter=maxiter)
    #@show stats
    sol = collect(sol) # to cpu
    params = reshape(sol,(p,N))
    chesse2mat.(eachcol(params),d)
end

"Container of the workflow results"
struct ddhResult{T<:AbstractFloat}
    # Basic parameters
    fdim::Int
    rdim::Int
    ssa::Bool
    # Vertex features
    u::Vector{T}
    div::Vector{T}
    rot::Vector{T}
    vgrass::Union{Vector{T},Nothing}
    spans::Vector{Matrix{T}}
    frames::Vector{Matrix{T}}
    planes::Vector{Matrix{T}}
    H::Vector{Matrix{T}}
    J::Vector{Matrix{T}}
    # Edge features
    w::Vector{T}
    egrass::Union{Vector{T},Nothing}
    rmaps::Union{Dict{Tuple{Int,Int}, Rmap}, Nothing}
    # Grad operators
    d0::SparseArrays.SparseMatrixCSC
    d0s::Union{SparseArrays.SparseMatrixCSC, Nothing}
end

function Base.show(io::IO, ::MIME"text/plain", ddh::ddhResult)
    fd, rd = size(ddh.spans[begin])
    ne, nv = size(ddh.d0)
    fd, rd, ssa = [getfield(ddh,x) for x in [:fdim,:rdim,:ssa]]
    ssa = ssa ? "on" : "off"
    print(io,":::ddHodgge::: values over graph {$nv, $ne}, dimensions: $fd -> $rd (ssa: $ssa)")
end

"""
    ddHodgeWorkflow(g, X, V;
        rdim=size(X,1)::Int, λ=0.1, ϵ=0.1,
        ssa=true, ssart=1, krytol=1e-32, useCUDA=false
    ) -> ddh::ddhResult

ddHodge standard workflow
"""
function ddHodgeWorkflow(
    g::SimpleGraph, X::AbstractMatrix, V::AbstractMatrix;
    rdim=size(X,1)::Int, λ=0.1::Float64, ϵ=λ::Float64,
    ssa=true::Bool, ssart=1::Int, krytol=1e-32::Float64, useCUDA=false::Bool
)
    fdim, N = size(X)
    # Easy timer
    elapsed(from) = @info "... done in $(round((time() - from),digits=3)) sec."
    @info "(1/4) Potential estimation"
    t0 = time()
    d0 = graphgrad(g) # gradient: C⁰ ⟶ C¹
    L0 = -d0'd0 # Δ₀
    w = calcwt(g,X,V) # edge weight ∈ C¹
    u, = linsolve(d0'd0,-d0'w,issymmetric=true) # potential ∈ C⁰
    elapsed(t0)
    @info "(2/4) Tangent space estimation"
    t0 = time()
    # Try to keep orientations of tangent spaces (local PCA) 
    spans = oriented_spans(g,X,rdim,ssart)
    concs = [sign(det(spans[src(e)]'spans[dst(e)])) for e in edges(g)]
    @info "Consistency of orientations: $(sum(concs .> 0)/ne(g))" 
    # Subspace alignment of local PCA spaces on nodes
    if ssa
        frames, d0s, rms = alignSpaces(g, spans; krytol=krytol, useCUDA=useCUDA)
        pas = hcat([rms[Tuple(e)].pa for e in edges(g)]...)
        elen = [norm(X[:,src(e)]-X[:,dst(e)]) for e in edges(g)]
        egrass = norm.(eachcol(pas)) # Grassmann distance
        vgrass = spdiagm(degree(g).^-1)*abs.(d0)'*(egrass./elen)
    else
        frames = @views spans
        d0s = nothing
        rms = nothing
        egrass = nothing
        vgrass = nothing
    end
    # Vertex-level Grassmann distance normalized by edge length
    elapsed(t0)
    @info "(3/4) Hessian matrix estimation"
    t0 = time()
    rhess = batchHess(g,X,V,L0,λ=λ,d=rdim,Ψ=frames,useCUDA=useCUDA)
    divr = -tr.(rhess) # divergence ∈ C⁰
    elapsed(t0)
    @info "(4/4) Jacobian matrix estimation"
    t0 = time()
    planes = [S[:,1:2] for S in spans]
    if ssa # also aligns PC1-2 planes
        planes, = alignSpaces(g, planes; krytol=krytol, useCUDA=false)
        # force to have same orientation in first 2 dim.
        fo = det(spans[ssart][:,1:2]'*planes[ssart])
        fo < 0 && reverse!.(planes,dims=2) # is costly operation?
    end
    oconcs = [sign(det(planes[src(e)]'planes[dst(e)])) for e in edges(g)]
    @info "Consistency of orientations (plane): $(sum(oconcs .> 0)/ne(g))" 
    Vdual = dualvels(V,planes,d=1:2) # dual velocities in the plane
    dhess = batchHess(g,X,Vdual,L0,λ=ϵ,d=2,Ψ=planes) # CUDA may not required
    rotr = tr.(dhess) # approx. rotation ∈ C₀ (dual)
    rS2W(r,S,W) = (W'W)\(W'*S)*[0 -r/2; r/2 0]*(S'*W) # rot in S to W
    #  rS2W(r,S,W) = pinv(W'W)*W'*S*[0 -r/2; r/2 0]*S'*W # rot in S to W
    # Reconstruct jacobian: J = D(F) = D(-grad) + D(rot) [+ D(harmonic)]
    rjacobs = [-rhess[i] + rS2W(rotr[i],planes[i],frames[i]) for i in 1:N]
    elapsed(t0)
    ddhResult{eltype(w)}(
        fdim, rdim, ssa,
        u, divr, rotr, vgrass, # vertex-level values
        spans, frames, planes, rhess, rjacobs, # vertex-level matrices
        w, egrass, rms, # edge-level features
        d0, d0s # operators
    )
end

"""
    schurselect(F::Schur, idx::Vector{Int}) -> (values, Z::Matrix)

Returns the basis of invariant subspace `Z`; the columns are the selected Schur vectors specified by `idx`.

!!! warning "Warning"
    The invariance is guaranteed if the conjugate pair of the complex eigenvalues were selected.

# Example
```jldoctest
julia> using LinearAlgebra: schur

julia> A = [0 1 0; -1 0 0; 0 0 1]
3×3 Matrix{Int64}:
  0  1  0
 -1  0  0
  0  0  1

julia> F = schur(A);

julia> i, j, = sortperm(abs.(imag.(F.values)),rev=true);

julia> vals, Z = schurselect(F,[i,j]) # x-y plane is the invariant of A 
(values = ComplexF64[0.0 + 1.0im, 0.0 - 1.0im], Z = [1.0 0.0; 0.0 1.0; 0.0 0.0])
```
"""
function schurselect(F::Schur, idx::Vector{Int})
    k = length(idx)
    T, Z, vals = ordschur(F, [i in idx for i in 1:length(F.values)])
    # first k vectors that span invariant subspace
    (values = vals[1:k], Z = Z[:,1:k])
end

"""
    basischange(A, V, W) -> B::Matrix

Change of basis: the matrix `A` in the basis `V` into the basis `W`

!!! warning "Warning"
    Use with caution when span(V) ≠ span(W).
"""
function basischange(A::AbstractMatrix, V::AbstractMatrix, W::AbstractMatrix)
    ((W'*W)\(W'*V))*A*((V'*V)\(V'*W))
end

"""
    basischange(A, V) -> B::Matrix

Change of basis: `A` in the basis `V` into the canonical basis, i.e., W = I
"""
function basischange(A::AbstractMatrix, V::AbstractMatrix)
    if size(A) != size(V)
        error("Error: change to cannonical basis is only allowed for the same size of square matrices")
    else
        V*A*((V'*V)\V')
    end
end

end # module ddHodge

