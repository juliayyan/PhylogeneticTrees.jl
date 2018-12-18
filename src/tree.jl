mutable struct BinaryTree
    depth::Int
    nnodes::Int
    edges::Vector{Tuple{Int,Int}}
    pathedges::Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}
    codes::Vector{Vector{Int}}
end 

function BinaryTree(depth::Int) 
    nnodes = 2^(depth+1)-1
    edges = Tuple{Int,Int}[]
    allcodes = Any[[-1],[0],[1]]
    codes = 0:1
    for d in 1:(depth-1) 
        codes = collect(IterTools.product(codes,0:1))
        push!(allcodes, codes)
    end
    allcodes = vcat([vec(c) for c in allcodes]...)
    function flatten(arr)
        rst = Any[]
        grep(v) =   for x in v
                    if isa(x, Tuple) 
                    grep(x) 
                    else push!(rst, x) end
                    end
        grep(arr)
        rst
    end
    allcodes = [reverse(flatten(code)) for code in allcodes]
    bt = BinaryTree(depth, nnodes, 
        edges,
        Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}(),
        allcodes)

    g = LightGraphs.Graph(nnodes)
    for d=0:(depth-1), n in getnodes(bt, d), c in getchildren(bt,n)
        push!(edges, (n,c))
        LightGraphs.add_edge!(g, (n,c))
    end
    fw = LightGraphs.floyd_warshall_shortest_paths(g)
    paths =  LightGraphs.enumerate_paths(fw) # [u][v]
    # (u,v) => list of edges on shortest path
    pathedges = Dict((u,v) => [getedge(bt, paths[u][v][i-1],paths[u][v][i]) 
                               for i in 2:length(paths[u][v])] 
                     for u in 1:nnodes for v in 1:nnodes)
    bt.pathedges = pathedges
    bt
end

validnode(bt::BinaryTree, node::Int) = ((node >= 1) && (node <= bt.nnodes))
validlayer(bt::BinaryTree, layer::Int) = ((layer >= 0) && (layer <= bt.depth))

function getlayer(bt::BinaryTree, node::Int)
    @assert validnode(bt, node)
    floor(Int, log(2,node))
end

function getnodes(bt::BinaryTree, layer::Int)
    #layer 0    1
    #layer 1    2:3
    @assert validlayer(bt, layer)
    (2^layer):(2^(layer+1)-1)
end 

getleaves(bt::BinaryTree) = getnodes(bt, bt.depth)

function getparent(bt::BinaryTree, node::Int)
    @assert validnode(bt, node)
    floor(Int, node/2)
end

function getchildren(bt::BinaryTree, node::Int)
    @assert validnode(bt, node)
    thislayer = getlayer(bt, node)
    thisid = findfirst(getnodes(bt, thislayer) .== node)
    nextnodes = getnodes(bt, thislayer+1)
    [nextnodes[2*thisid-1], nextnodes[2*thisid]]
end

function getdescendants(bt::BinaryTree, node::Int)
    @assert validnode(bt, node)
    thislayer = getlayer(bt, node)
    descendants = Int[]
    if thislayer < bt.depth 
        descendants = getchildren(bt, node)
        for d in descendants
            dlayer = getlayer(bt, d)
            descendants = vcat(descendants, getdescendants(bt, d))
        end
    end
    return descendants
end

function getsubtreeleaves(bt::BinaryTree, node::Int)
    return intersect(vcat(node,getdescendants(bt, node)), getleaves(bt))
end

function getedge(bt::BinaryTree, u::Int, v::Int)
    @assert validnode(bt, u)
    @assert validnode(bt, v)
    return (min(u,v),max(u,v))
end

function Base.show(io::IO, bt::BinaryTree; offset::String="")
    println(io, offset, string(typeof(bt)))
    println(io, offset, "    * depth: $(bt.depth)  (nnodes: $(bt.nnodes))")
end
