mutable struct CodedTreeProblem
    pd::PopulationData
    bt::BinaryTree
    outgroupnode::Int
    model::JuMP.Model
    assign::JuMP.JuMPArray{JuMP.Variable}
    countedge::JuMP.JuMPDict{JuMP.Variable}
    weight::JuMP.JuMPArray{JuMP.Variable}
    f3formula
    f3err::JuMP.JuMPDict{JuMP.Variable}
end 

function CodedTreeProblem(
    pd::PopulationData, 
    bt::BinaryTree;
    solver = Gurobi.GurobiSolver(),
    binaryencoding::Bool = false)
    
    const npop   = pd.npop
    const edges  = bt.edges
    const leaves = getleaves(bt)
    const outgroupnode = leaves[1]

    tree = JuMP.Model(solver=solver)
    if binaryencoding 
        JuMP.@variable(tree, assign[1:npop,getleaves(bt)] >= 0)
        JuMP.@variable(tree, codeselect[1:pd.npop,1:bt.depth], Bin)
        JuMP.@constraint(tree, 
            [a=1:pd.npop,m=1:bt.depth], 
            codeselect[a,m] == sum(bt.codes[u][m]*assign[a,u] for u in leaves))

    else 
        JuMP.@variable(tree, assign[1:npop,getleaves(bt)], Bin)
    end
    JuMP.@variable(tree, weight[edges] >= 0)
    JuMP.@variable(tree, weightaux[a=1:npop,b=a:npop,edges] >= 0)
    JuMP.@variable(tree, countedge[a=1:npop,b=a:npop,edges] >= 0)
    JuMP.@expression(tree, 
        f3formula[a=1:npop,b=a:npop],
        sum(weightaux[a,b,edg] for edg in edges))
    JuMP.@variable(tree, f3err[a=1:npop,b=a:npop])

    validtreeconstraints(pd, bt, tree, assign, outgroupnode, nlevels)
    countedgeconstraints(pd, bt, tree, assign, countedge, outgroupnode)
    errorconstraints(pd, bt, tree, weight, weightaux, countedge, f3formula, f3err)

    JuMP.@objective(tree, Min, 
        sum(pd.cov[a1,b1,a2,b2]*f3err[a1,b1]*f3err[a2,b2] 
            for a1 in 1:npop, a2 in 1:npop, b1 in a1:npop, b2 in a2:npop))

    CodedTreeProblem(
        pd, bt, outgroupnode, 
        tree, 
        assign, countedge, 
        weight, 
        f3formula, f3err)
end

function countedgeconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::JuMP.JuMPArray{JuMP.Variable},
    countedge::JuMP.JuMPDict{JuMP.Variable},
    outgroupnode::Int)
    
    const npop = pd.npop
    
    for (u,v) in bt.edges 
        if in((u,v), bt.pathedges[outgroupnode,1])
            JuMP.@constraint(tree, 
                [a=1:npop, b=a:npop],
                countedge[a,b,(u,v)] <= 1-sum(assign[a,n] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop],
                countedge[a,b,(u,v)] <= 1-sum(assign[b,n] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop],
                countedge[a,b,(u,v)] >= 1 - sum(assign[a,n] + assign[b,n] for n in getsubtreeleaves(bt,v)))
        else
            JuMP.@constraint(tree, 
                [a=1:npop, b=a:npop],
                countedge[a,b,(u,v)] <= sum(assign[a,n] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop],
                countedge[a,b,(u,v)] <= sum(assign[b,n] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop],
                countedge[a,b,(u,v)] >= sum(assign[a,n] + assign[b,n] for n in getsubtreeleaves(bt,v)) - 1)
        end
    end

end

function errorconstraints(pd::PopulationData, 
    bt::BinaryTree, 
    tree::JuMP.Model, 
    weight::JuMP.JuMPArray{JuMP.Variable}, 
    weightaux::JuMP.JuMPDict{JuMP.Variable}, 
    countedge::JuMP.JuMPDict{JuMP.Variable}, 
    f3formula,#::JuMP.JuMPDict{JuMP.Variable}, 
    f3err::JuMP.JuMPDict{JuMP.Variable})

    const npop = pd.npop 
    const bigm = maximum(pd.f3)*2

    # set weight bilinear terms
    JuMP.@constraint(tree, 
        [a=1:npop,b=a:npop,edg=bt.edges],
        weightaux[a,b,edg] <= bigm*countedge[a,b,edg])
    JuMP.@constraint(tree, 
        [a=1:npop,b=a:npop,edg=bt.edges],
        weightaux[a,b,edg] <= weight[edg])
    JuMP.@constraint(tree, 
        [a=1:npop,b=a:npop,edg=bt.edges],
        weightaux[a,b,edg] >= weight[edg] + bigm*countedge[a,b,edg] - bigm)

    # set error terms
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] >= pd.f3[a,b] - f3formula[a,b])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] >= f3formula[a,b] - pd.f3[a,b])

end

# warning: lots of magic constants here
function printtree(tp::CodedTreeProblem)

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
            if d == depth
                nodeassign = round.(JuMP.getvalue(tp.assign[:,n]))
            else 
                nodeassign = 0
            end
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

function Base.show(io::IO, tp::CodedTreeProblem; offset::String="")
    println(io, offset, string(typeof(tp)))
end
