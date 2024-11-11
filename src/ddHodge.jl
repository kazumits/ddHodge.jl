module ddHodge

export kNNGraph, ddHodgeWorkflow

using LinearAlgebra, SparseArrays, Graphs
using KrylovKit: linsolve, eigsolve
using NearestNeighbors: KDTree, knn


module DEC
    function delaunayGraph end
    function stardiags end
end

module GPU
    function cuspm end
    function cuvec end
end

"""
    graphgrad(g)

Construct grad operator from graph `g`

The direction is implicitly determined by node order, i.e., `i -> j` if `i < j`.

# Example
```jldoctest
julia> using ddHodge.Graphs

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
function calcwt(g::SimpleGraph{Int},pts::AbstractMatrix,vel::AbstractMatrix)
    map(edges(g)) do e
        i, j = src(e), dst(e)
        dx = pts[:,j] - pts[:,i] # i to j
        (vel[:,i] + vel[:,j])'*dx/2.0
    end
end

"""
    kNNGraph(X, k)
Graph construction using k-NN
"""
function kNNGraph(X::AbstractMatrix,k::Int)
    N = size(X,2)
    kdtree = KDTree(X)
    idxs, dists = knn(kdtree, X, k+1, true)
    g = SimpleGraph(N)
    for (i,v) in enumerate(idxs)
        for j in 2:(k+1)
            # disconnect if zero distance
            if dists[i][j] > 0.0
                add_edge!(g,v[1],v[j])
            end
        end
    end
    #rem_vertices!(g,findall(degree(g).==0))
    #@assert is_connected(g) # stop when g is not connected
    g, kdtree, idxs, dists
end
 
# Coefficients of directional derivatives w.r.t. Hessian
function chesse(v)
    n = length(v)
    c = fill(0.0::Float64,Int(n*(n-1)/2))
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
function chesse2mat(x,d)
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
function dualvels(V, spans; d=[1,2])
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
function getspan(g::SimpleGraph{Int},pts::AbstractMatrix,i::Int,d::Int; ref=nothing)
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
function oriented_spans(g::SimpleGraph{Int},pts::AbstractMatrix,d::Int,root::Int)
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

# Batch Hessian matrix estimation
function batchHessCU(
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
    p = size(Ms[1],2)
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


"""
    ddHodgeWorkflow(g, X, V;
        ldim=size(X,1)::Int, λ=0.1, ϵ=0.1, ssa=true, ssart=1, krytol=1e-32, useCUDA=false)

ddHodge standard workflow
"""
function ddHodgeWorkflow(
    g::SimpleGraph, X::AbstractMatrix, V::AbstractMatrix;
    ldim=size(X,1)::Int, λ=0.1, ϵ=0.1, ssa=true, ssart=1, krytol=1e-32, useCUDA=false
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
    #spans = fdim > ldim ? oriented_spans(g,X,ldim,ssart) : fill(Matrix(1.0I,fdim,fdim),N)
    spans = oriented_spans(g,X,ldim,ssart)
    concs = [sign(det(spans[src(e)]'spans[dst(e)])) for e in edges(g)]
    @info "Consistency of orientations: $(sum(concs .> 0)/ne(g))" 
    # Restriction maps and sheaf on graph
    rms = Dict([Tuple(e) => connectspans(spans[src(e)],spans[dst(e)]) for e in edges(g)])
    d0s = d0sheaf(g,rms) # Sheaf gradient
    if ssa
        cuspm = useCUDA ? GPU.cuspm : identity
        cuvec = useCUDA ? GPU.cuvec : identity
        L0s = cuspm(d0s'd0s)
        x0 = cuvec(rand(N*ldim))
        #if usenormalL # faster calculation for Hessian?
            #Dhi = spdiagm(repeat(degree(g).^(-1/2),inner=ldim))
            #L0s = Dhi*L0s*Dhi
        #end
        vals, vecs, = eigsolve(L0s,x0,ldim,:SR,issymmetric=true,tol=krytol)
        vecs = collect.(vecs) # to cpu
        vec2dir(x) = normalize.(eachcol(reshape(x,(ldim,N))))
        # Basis matrix of each local dimension
        frames = [hcat(spans.*vec2dir(v)...) for v in vecs[1:ldim]]
        # Basis matrix of each vertex
        Ψ = @views [hcat([frames[i][:,j] for i in 1:ldim]...) for j in 1:N]
    else
        Ψ = @views spans
    end
    # Vertex-level Grassmann distance normalized by edge length
    pas = hcat([rms[Tuple(e)].pa for e in edges(g)]...)
    elen = [norm(X[:,src(e)]-X[:,dst(e)]) for e in edges(g)]
    vgrass = spdiagm(degree(g).^-1)*abs.(d0)'*(norm.(eachcol(pas))./elen)
    elapsed(t0)
    @info "(3/4) Hessian matrix estimation"
    t0 = time()
    rhess = batchHessCU(g,X,V,L0,λ=λ,d=ldim,Ψ=Ψ,useCUDA=useCUDA)
    divr = -tr.(rhess) # divergence ∈ C⁰
    elapsed(t0)
    @info "(4/4) Jacobian matrix estimation"
    t0 = time()
    if ssa
        rplanes = [S[:,1:2] for S in spans]
        ofix(cur,ref) = det(cur'*ref) < 0 ? [cur[:,1:end-1] -cur[:,end]] : cur
        oplanes = propagate((pa,ch,ref) -> ofix(rplanes[ch],ref), g, rplanes[ssart], root=ssart)
    else
        oplanes = oriented_spans(g,X,2,ssart)
    end
    oconcs = [sign(det(oplanes[src(e)]'oplanes[dst(e)])) for e in edges(g)]
    @info "Consistency of orientations (plane): $(sum(oconcs .> 0)/ne(g))" 
    Vdual = dualvels(V,oplanes,d=1:2) # dual velocities in PC1-2 plane
    dhess = batchHessCU(g,X,Vdual,L0,λ=ϵ,d=2,Ψ=oplanes) # CUDA may not required
    rotr = tr.(dhess) # approx. rotation ∈ C₀ (dual)
    R2Ψ(r,S,W) = pinv(W'W)*W'*S*[0 -r/2; r/2 0]*S'*W # A in S to A in Ψ
    # reconstruction of jacobian: J = D(F) = D(-grad) + D(rot) + D(harmonic)
    rjacobs = [-rhess[i] + R2Ψ(rotr[i],oplanes[i],Ψ[i]) for i in 1:N]
    elapsed(t0)
    (
        ldim=ldim,
        w=w, u=u, div=divr, rot=rotr, vgrass=vgrass,
        d0=d0, d0s=d0s, spans=spans, oplanes=oplanes,
        rmaps=rms, frames=Ψ,H=rhess, J=rjacobs
    )
end

"""
    schursel(F::Schur, idx)

Get selected invariant subspace calculated by `schur`. 
"""
function schursel(F::Schur,idx::Vector{Int})
    k = length(idx)
    T, Z, vals = ordschur(F,[i in idx for i in 1:length(F.values)])
    vals[1:k], Z[:,1:k] # first k vectors span invariant sub.
end

"""
    basistrans(X, V, W)

Basis trans. of matrix `X` in `V` into `W`

Caution: dim(V ⋂ W) should not be 0.
"""
basistrans(X::AbstractMatrix,V::AbstractMatrix,W::AbstractMatrix) = ((W'*W)\(W'*V))*X*((V'*V)\(V'*W))

end # module ddHodge

