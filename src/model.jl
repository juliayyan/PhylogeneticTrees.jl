mutable struct TreeProblem
    pd::PopulationData
    bt::BinaryTree
    outgroupnode::Int
    nlevels::Int
    model::JuMP.Model
    assign::JuMP.JuMPArray{JuMP.Variable}
    assign2::JuMP.JuMPDict{JuMP.Variable}
    weight::JuMP.JuMPArray{JuMP.Variable}
    f3formula::JuMP.JuMPDict{JuMP.Variable}
    f3err::JuMP.JuMPDict{JuMP.Variable}
end 

function TreeProblem(
    pd::PopulationData, 
    bt::BinaryTree;
    outgroupnode::Int = 2^bt.depth,
    binaryencoding::Bool = false,
    nlevels::Int = 1,
    solver = Gurobi.GurobiSolver())
    
    @assert in(outgroupnode, getnodes(bt,bt.depth))
    @assert !binaryencoding || (nlevels == 1)

    const npop   = pd.npop
    const edges  = bt.edges
    const othernodes = setdiff(getnodes(bt,bt.depth), outgroupnode)
    const levels = 1:nlevels

    tree = JuMP.Model(solver=solver)
    if !binaryencoding
        JuMP.@variable(tree, assign[1:npop,getnodes(bt,bt.depth),levels], Bin)
    else 
        JuMP.@variable(tree, assign[1:npop,getnodes(bt,bt.depth),levels] >= 0)
    end
    JuMP.@variable(tree, assign2[a=1:npop,b=a:npop,othernodes,othernodes,levels,levels])
    JuMP.@variable(tree, weight[edges] >= 0)
    JuMP.@variable(tree, f3formula[a=1:npop,b=a:npop,u=othernodes,v=othernodes,levels,levels] >= 0)
    JuMP.@variable(tree, f3err[a=1:npop,b=a:npop])
 
    validtreeconstraints(pd, bt, tree, assign, outgroupnode, nlevels)
    logicalconstraints(pd, bt, tree, assign, assign2, outgroupnode, nlevels)
    errorconstraints(pd, bt, tree, assign2, weight, f3formula, f3err, outgroupnode, nlevels)
    binaryencoding && binaryencodingconstraints(pd, bt, tree, assign)

    JuMP.@objective(tree, Min, 
        sum(pd.cov[a1,b1,a2,b2]*f3err[a1,b1]*f3err[a2,b2] 
            for a1 in 1:npop, a2 in 1:npop, b1 in a1:npop, b2 in a2:npop))

    TreeProblem(pd, bt, outgroupnode, nlevels, tree, assign, assign2, weight, f3formula, f3err)
end

function validtreeconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::JuMP.JuMPArray{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int)

    const npop = pd.npop
    const levels = 1:nlevels
    const leaves = getnodes(bt,bt.depth)
    
    # one-to-one assignment
    JuMP.@constraint(tree, 
        [a=1:npop], 
        sum(assign[a,leaves,levels]) == nlevels)
    JuMP.@constraint(tree, 
        [n=getnodes(bt,bt.depth)], 
        sum(assign[1:npop,n,1]) <= 1)
    JuMP.@constraint(tree, 
        [a=1:npop,n=leaves,l=2:nlevels], 
        assign[a,n,l] <= assign[a,n,l-1])
    # assign outgroup to a node
    JuMP.@constraint(tree, 
        [l=levels], 
        assign[pd.outgroup,outgroupnode,l] == 1)

end

function logicalconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::JuMP.JuMPArray{JuMP.Variable},
    assign2::JuMP.JuMPDict{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int)

    const npop   = pd.npop
    const othernodes = setdiff(getnodes(bt,bt.depth), outgroupnode)
    const levels = 1:nlevels
    
    JuMP.@constraint(tree, [a=1:npop,b=a:npop,u=othernodes,v=othernodes,l=levels,m=levels],
        assign2[a,b,u,v,l,m] <= assign[a,u,l])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop,u=othernodes,v=othernodes,l=levels,m=levels],
        assign2[a,b,u,v,l,m] <= assign[b,v,m])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop,u=othernodes,v=othernodes,l=levels,m=levels],
        assign2[a,b,u,v,l,m] >= assign[a,u,l] + assign[b,v,m] - 1)

end 

function errorconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign2::JuMP.JuMPDict{JuMP.Variable},
    weight::JuMP.JuMPArray{JuMP.Variable},
    f3formula::JuMP.JuMPDict{JuMP.Variable},
    f3err::JuMP.JuMPDict{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int)

    const npop   = pd.npop
    const othernodes = setdiff(getnodes(bt,bt.depth), outgroupnode)
    const bigm = maximum(pd.f3)*2
    const levels = 1:nlevels

    JuMP.@expression(tree,
        f3pathsum[u=othernodes,v=othernodes],
        sum(weight[edg] for edg in 
            intersect(bt.pathedges[u,outgroupnode], 
                      bt.pathedges[v,outgroupnode])))
    JuMP.@constraint(tree,
        [a=1:npop,b=a:npop,u=othernodes,v=othernodes,l=levels,m=levels],
        f3formula[a,b,u,v,l,m] <= bigm*assign2[a,b,u,v,l,m])
    JuMP.@constraint(tree,
        [a=1:npop,b=a:npop,u=othernodes,v=othernodes,l=levels,m=levels],
        f3formula[a,b,u,v,l,m] >= f3pathsum[u,v] + bigm*assign2[a,b,u,v,l,m] - bigm)
    JuMP.@constraint(tree,
        [a=1:npop,b=a:npop,u=othernodes,v=othernodes,l=levels,m=levels],
        f3formula[a,b,u,v,l,m] <= f3pathsum[u,v])
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] >= pd.f3[a,b] - 
                      sum(f3formula[a,b,u,v,l,m] 
                          for u=othernodes,v=othernodes,l=levels,m=levels)/nlevels^2)
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] >= sum(f3formula[a,b,u,v,l,m] 
                          for u=othernodes,v=othernodes,l=levels,m=levels)/nlevels^2 - 
                      pd.f3[a,b])

end 

# Juan Pablo trickery to improve branch-and-bound performance
function binaryencodingconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model,
    assign::JuMP.JuMPArray{JuMP.Variable})
    const leaves = getnodes(bt,bt.depth)
    codes = 0:1
    dim = bt.depth
    for d in 1:(dim-1) 
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
    JuMP.@variable(tree, codeselect[1:pd.npop,1:dim], Bin)
    JuMP.@constraint(tree, 
        [a=1:pd.npop,m=1:dim], 
        codeselect[a,m] == sum(codes[i][m]*assign[a,leaves[i],1] for i in 1:length(leaves)))
end

function warmstartunmixed(tp::TreeProblem; 
    timelimit::Int = 30,
    solver = Gurobi.GurobiSolver(OutputFlag = 0, TimeLimit = timelimit))
    @assert tp.nlevels > 1
    tp0 = TreeProblem(tp.pd, tp.bt, solver = solver)
    breaksymmetries(tp0)
    JuMP.solve(tp0.model)
    solution0 = JuMP.getvalue(tp0.assign);
    for a in 1:tp.pd.npop, u in getnodes(tp.bt, tp.bt.depth), l in 1:tp.nlevels 
        JuMP.setvalue(tp.assign[a,u,l], JuMP.getvalue(tp0.assign[a,u,1]))
    end
end

function printnodes(tp::TreeProblem)
    for u in getnodes(tp.bt, tp.bt.depth), a in 1:tp.pd.npop
        level = round(sum(JuMP.getvalue(tp.assign[a,u,:])))/tp.nlevels
        level > 0 && println(tp.pd.pops[a][1:2], "\t", u, "\t", level)
    end
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
            if d == depth
                nodeassign = round.(JuMP.getvalue(tp.assign[:,n,1]))
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

function Base.show(io::IO, tp::TreeProblem; offset::String="")
    println(io, offset, string(typeof(tp)))
end
