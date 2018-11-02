mutable struct BinaryTree
    depth::Int
    nnodes::Int
    edges::Vector{Tuple{Int,Int}}
    pathedges::Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}
end 

function BinaryTree(depth::Int) 
    nnodes = 2^(depth+1)-1
    edges = Tuple{Int,Int}[]
    bt = BinaryTree(depth, nnodes, 
        edges,
        Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}())

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

function getedge(bt::BinaryTree, u::Int, v::Int)
    @assert validnode(bt, u)
    @assert validnode(bt, v)
    return (min(u,v),max(u,v))
end

function Base.show(io::IO, bt::BinaryTree; offset::String="")
    println(io, offset, string(typeof(bt)))
    println(io, offset, "    * depth: $(bt.depth)  (nnodes: $(bt.nnodes))")
end
