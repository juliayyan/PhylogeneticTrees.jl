function validtreeconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::JuMP.JuMPArray{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int)

    const npop = pd.npop
    const levels = 1:nlevels
    const leaves = getleaves(bt)
    
    # one-to-one assignment
    JuMP.@constraint(tree, 
        [a=1:npop], 
        sum(assign[a,leaves,levels]) == nlevels)
    JuMP.@constraint(tree, 
        [n=getleaves(bt)], 
        sum(assign[1:npop,n,1]) <= 1)
    JuMP.@constraint(tree, 
        [a=1:npop,n=leaves,l=2:nlevels], 
        assign[a,n,l] <= assign[a,n,l-1])
    # assign outgroup to a node
    JuMP.@constraint(tree, 
        [l=levels], 
        assign[pd.outgroup,outgroupnode,l] == 1)

end

# Juan Pablo trickery to improve branch-and-bound performance
function binaryencodingconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model,
    assign::JuMP.JuMPArray{JuMP.Variable})

    JuMP.@variable(tree, codeselect[1:pd.npop,1:bt.depth], Bin)
    JuMP.@constraint(tree, 
        [a=1:pd.npop,m=1:bt.depth], 
        codeselect[a,m] == sum(bt.codes[u][m]*assign[a,u,1] for u in getleaves(bt)))

end

function warmstartunmixed(tp::Union{NodeTreeProblem,TreeProblem}; 
    timelimit::Int = 30,
    solver = Gurobi.GurobiSolver(OutputFlag = 0, TimeLimit = timelimit))
    @assert tp.nlevels > 1
    tp0 = TreeProblem(tp.pd, tp.bt, solver = solver)
    JuMP.solve(tp0.model)
    solution0 = JuMP.getvalue(tp0.assign);
    for a in 1:tp.pd.npop, u in getnodes(tp.bt, tp.bt.depth), l in 1:tp.nlevels 
        JuMP.setvalue(tp.assign[a,u,l], JuMP.getvalue(tp0.assign[a,u,1]))
    end
end

# fixes population a to be un-admixed
function unmix(tp::Union{NodeTreeProblem,TreeProblem}, a::Int)
    JuMP.@constraint(tp.model, sum(tp.assign[a,:,1]) == 1)
end

# fixes population a to node u
function fix(tp::Union{NodeTreeProblem,TreeProblem}, a::Int, u::Int)
    JuMP.@constraint(tp.model, 
        [l=1:tp.nlevels],
        tp.assign[a,u,l] == 1)
end

function removesolution(
    tp::Union{NodeTreeProblem,TreeProblem},
    solution::JuMP.JuMPArray{Float64})
    expr = 0
    for a in 1:tp.pd.npop, n in getnodes(tp.bt, tp.bt.depth), l in 1:tp.nlevels
        if solution[a,n,l] > 0
            expr += tp.assign[a,n,l]
        else 
            expr += (1-tp.assign[a,n,l])
        end
    end
    JuMP.@constraint(tp.model, expr <= JuMP.getvalue(expr)-1)
end
