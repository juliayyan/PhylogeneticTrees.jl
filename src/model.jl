mutable struct TreeProblem
    pd::PopulationData
    bt::BinaryTree
    outgroupnode::Int
    model::JuMP.Model
    assign::Matrix{JuMP.Variable}
    assign2::JuMP.JuMPDict{JuMP.Variable,4}
    weight::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{Tuple{Int64,Int64},1}}}
    f3formula::JuMP.JuMPDict{JuMP.Variable,4}
    f3err::JuMP.JuMPDict{JuMP.Variable,2}
end 

function TreeProblem(
    pd::PopulationData, 
    bt::BinaryTree,
    outgroupnode::Int;
    solver = Gurobi.GurobiSolver())

    const npop   = pd.npop
    const nnodes = bt.nnodes
    const edges  = bt.edges
    const othernodes = setdiff(2:nnodes, outgroupnode)

    tree = JuMP.Model(solver=solver)
    JuMP.@variable(tree, assign[1:npop,1:nnodes], Bin) # switch to >= 0 if binary encoding
    JuMP.@variable(tree, assign2[a=1:npop,b=a:npop,othernodes,othernodes])
    JuMP.@variable(tree, weight[edges] >= 0)
    JuMP.@variable(tree, f3formula[a=1:npop,b=a:npop,u=othernodes,v=othernodes] >= 0)
    JuMP.@variable(tree, f3err[a=1:npop,b=a:npop])
 
    validtreeconstraints(pd, bt, tree, assign, outgroupnode)
    logicalconstraints(pd, bt, tree, assign, assign2, outgroupnode)
    errorconstraints(pd, bt, tree, assign2, weight, f3formula, f3err, outgroupnode)
    #binaryencodingconstraints(pd, bt, tree, assign)

    JuMP.@objective(tree, Min, 
        sum(pd.cov[a1,b1,a2,b2]*f3err[a1,b1]*f3err[a2,b2] 
            for a1 in 1:npop, a2 in 1:npop, b1 in a1:npop, b2 in a2:npop))

    TreeProblem(pd, bt, outgroupnode, tree, assign, assign2, weight, f3formula, f3err)
end

function validtreeconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::Matrix{JuMP.Variable},
    outgroupnode::Int)

    const npop   = pd.npop
    const nnodes = bt.nnodes
    
    # nothing at root
    JuMP.@constraint(tree, sum(assign[1:npop,1]) == 0)
    # one-to-one assignment
    JuMP.@constraint(tree, [a=1:npop],   sum(assign[a,1:nnodes]) == 1)
    JuMP.@constraint(tree, [n=1:nnodes], sum(assign[1:npop,n])   <= 1)
    # populations are only assigned to leaf nodes
    JuMP.@constraint(tree, [d=0:(bt.depth-1),
                            n=getnodes(bt,d),
                            c=getdescendants(bt,n)], 
        sum(assign[1:npop,c]) <= 1-sum(assign[1:npop,n]))
    # assign outgroup to a node
    JuMP.@constraint(tree, assign[pd.outgroup,outgroupnode] == 1)

end

function logicalconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::Matrix{JuMP.Variable},
    assign2::JuMP.JuMPDict{JuMP.Variable,4},
    outgroupnode::Int
    )

    const npop   = pd.npop
    const nnodes = bt.nnodes
    const othernodes = setdiff(2:nnodes, outgroupnode)
    
    JuMP.@constraint(tree, [a=1:npop,b=a:npop,u=othernodes,v=othernodes],
        assign2[a,b,u,v] <= assign[a,u])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop,u=othernodes,v=othernodes],
        assign2[a,b,u,v] <= assign[b,v])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop,u=othernodes,v=othernodes],
        assign2[a,b,u,v] >= assign[a,u] + assign[b,v] - 1)

end 

function errorconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign2::JuMP.JuMPDict{JuMP.Variable,4},
    weight::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{Tuple{Int64,Int64},1}}},
    f3formula::JuMP.JuMPDict{JuMP.Variable,4},
    f3err::JuMP.JuMPDict{JuMP.Variable,2},
    outgroupnode::Int)

    const npop   = pd.npop
    const nnodes = bt.nnodes
    const othernodes = setdiff(2:nnodes, outgroupnode)
    const bigm = maximum(pd.f3)*2

    JuMP.@expression(tree,
        f3pathsum[u=othernodes,v=othernodes],
        sum(weight[edg] for edg in 
            intersect(bt.pathedges[u,outgroupnode], 
                      bt.pathedges[v,outgroupnode])))
    JuMP.@constraint(tree,
        [a=1:npop,b=a:npop,u=othernodes,v=othernodes],
        f3formula[a,b,u,v] <= bigm*assign2[a,b,u,v])
    JuMP.@constraint(tree,
        [a=1:npop,b=a:npop,u=othernodes,v=othernodes],
        f3formula[a,b,u,v] >= f3pathsum[u,v] + bigm*assign2[a,b,u,v] - bigm)
    JuMP.@constraint(tree,
        [a=1:npop,b=a:npop,u=othernodes,v=othernodes],
        f3formula[a,b,u,v] <= f3pathsum[u,v])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] >= pd.f3[a,b] - sum(f3formula[a,b,u,v] for u=othernodes,v=othernodes))
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] >= sum(f3formula[a,b,u,v] for u=othernodes,v=othernodes) - pd.f3[a,b])

end 

# Juan Pablo trickery to improve branch-and-bound performance
function binaryencodingconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model,
    assign::Matrix{JuMP.Variable})
    # speedup attempts
    # symmetry breaking constraints don't seem to be helpful
    # preventing chains of unweighted edges doesn't seem to be helpful
    codes = 0:1
    dim = bt.depth
    for d in 1:dim 
        codes = collect(Iterators.product(codes,0:1))
    end
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
    codes = [flatten(code) for code in codes]
    JuMP.@variable(tree, codeselect[1:pd.npop,1:(dim+1)], Bin)
    JuMP.@constraint(tree, 
        [a=1:pd.npop,m=1:(dim+1)], 
        codeselect[a,m] == sum(codes[n][m]*assign[a,n] for n in 1:bt.nnodes))
end

# warning: lots of magic constants here
function printtree(tp::TreeProblem)

    const pd = tp.pd
    const bt = tp.bt
    const depth = bt.depth

    spaces = Dict(zip(0:depth, [3 for i in 0:depth]))
    for d=(depth-1):-1:0
        spaces[d] = 2*spaces[d+1]+1
    end
    for d=0:depth
        nodes = getnodes(bt,d)
        for n in nodes 
            print(" "^spaces[d])
            nodeassign = round.(JuMP.getvalue(tp.assign[:,n]))
            if sum(nodeassign) < 1 
                str = "()"
            else 
                str = pd.pops[findfirst(nodeassign)][1:2]
            end
            @printf("%2s", str)
            print(" "^spaces[d])
        end
        println()
        if d < depth 
            for n in nodes 
                print(" "^(spaces[d]-2), " /"," "^2,"\\ "," "^(spaces[d]-2))
            end
            println()
            for n in nodes 
                weight1 = JuMP.getvalue(tp.weight[(n,getchildren(bt,n)[1])])
                weight2 = JuMP.getvalue(tp.weight[(n,getchildren(bt,n)[2])])
                print(" "^(spaces[d]-3))
                if round(weight1) > 0
                    @printf("%3.0f", weight1)
                else 
                    @printf("%3s", " ")
                end
                print(" "^2)
                if round(weight2) > 0
                    @printf("%3.0f", weight2)
                else 
                    @printf("%3s", " ")
                end
                print(" "^(spaces[d]-3))
            end
            println()
            for n in nodes 
                print(" "^(spaces[d]-2),"/ "," "^2," \\"," "^(spaces[d]-2))
            end
            println()
        end  
    end

end

function Base.show(io::IO, tp::TreeProblem; offset::String="")
    println(io, offset, string(typeof(tp)))
end
