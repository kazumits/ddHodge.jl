module DECExt

export delaunayGraph, stardiags

using Triangulate, Graphs, SparseArrays, LinearAlgebra
using ddHodge.DEC

"""
    DEC.delaunayGraph(pts; option="Qev")

Construct graph using Delaunay triangulation and the curl operator
"""
function DEC.delaunayGraph(pts::AbstractMatrix;option="Qev")
    nip = size(pts,2)
    triio = TriangulateIO()
    triio.pointlist = pts
    tout, vout = triangulate(option,triio)
    # primal to dual edges
    duale = Dict(
        zip(
            Tuple.(sort.(eachcol(tout.edgelist))),
            Tuple.(sort.(eachcol(vout.edgelist)))
        )
    )
    filter!(p -> p.second[1] > 0, duale)
    etri = tout.edgelist
    vtri = sort(tout.trianglelist,dims=1)
    nop = maximum(vtri)
    h = Graph(nop) # may contains steiner points
    for i in 1:size(etri,2)
        add_edge!(h,etri[1,i],etri[2,i])
    end
    nip < nop && rem_vertices!(h,collect(nip+1:nop))
    # remove if the triangle contains steiner points
    vtri = vtri[:,vec(sum(vtri .<= nip,dims=1) .== 3)]
    ntri = size(vtri,2)
    edic = Dict(zip(edges(h),1:ne(h)))
    eidx = (i,j) -> edic[Edge(i,j)] # to get edge idxs
    itri = mapslices(x -> eidx.(x[[1,2,1]],x[[2,3,3]]),vtri,dims=1)
    # [1,1,-1] since 1 -> 3 is the positive direction in graph
    curl = sparse(repeat(1:ntri,inner=3),vec(itri),repeat([1,1,-1],ntri))
    return(graph=h,triangles=vtri,curl=curl,voronoi=vout,duale=duale)
end

"""
    DEC.stardiags(g, pts, tri)

Returns diagonal elements of Hodge star of triangle mesh
"""
function DEC.stardiags(g::SimpleGraph,pts::AbstractMatrix,tri::Matrix{Int32})
    # cotangents of triangles
    cotans = similar(tri,Float64)
    for t in 1:size(tri,2)
        x = @views pts[:,tri[:,t]]
        for (i,j,k) in [(1,2,3),(2,3,1),(3,1,2)]
            v = normalize(x[:,j] - x[:,i])
            w = normalize(x[:,k] - x[:,i])
            cosθ  = clamp(v'w,-1.0,1.0)
            cotans[i,t] = cosθ /sqrt(1-cosθ ^2)
        end
    end

    # dual 0-form
    vdual = fill(0.0,nv(g))
    for t in 1:size(tri,2)
        v = @views tri[:,t]
        x = @views pts[:,v]
        for (i,j,k) in [(1,2,3),(2,3,1),(3,1,2)]
            Lij = norm(x[:,j] - x[:,i])
            Lik = norm(x[:,k] - x[:,i])
            area = (Lik^2*cotans[j,t] + Lij^2*cotans[k,t])/8
            vdual[v[i]] += area
        end
    end

    # dual 1-form
    d1scale = Dict{Tuple{Int,Int},Float64}()
    for t in 1:size(tri,2)
        v = tri[:,t]
        for (i,j,k) in [(1,2,3),(2,3,1),(3,1,2)]
            s, d = (v[i] < v[j]) ? (v[i], v[j]) : (v[j], v[i])
            d1scale[(s,d)] = cotans[k,t]/2 + get(d1scale,(s,d),0.0)
        end
    end
    edual = [get(d1scale,(src(e),dst(e)),0.0) for e in edges(g)]

    return vdual, edual
end

end # module DECExt
