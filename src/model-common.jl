function validtreeconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::JuMP.JuMPArray{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int,
    nmixtures::Int)

    npop = pd.npop
    levels = 1:nlevels
    leaves = getleaves(bt)
    
    # each population gets assigned levels
    JuMP.@constraint(tree,
        [a=1:npop,l=levels],
        sum(assign[a,leaves,l]) == 1)
    # each node gets at most one population
    if nlevels > 1
        JuMP.@variable(tree, assign1[1:npop,leaves] >= 0)
        JuMP.@constraint(tree,
            [a=1:npop,n=leaves,l=levels],
            assign1[a,n] >= assign[a,n,l])
        JuMP.@constraint(tree,
            [n=leaves],
            sum(assign1[a,n] for a in 1:npop) <= 1)
        JuMP.@constraint(tree,
            sum(assign1) <= nmixtures)
    else 
        JuMP.@constraint(tree,
            [n=leaves],
            sum(assign[1:npop,n,1]) <= 1)
    end
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
        JuMP.setvalue(tp.model[:assign1][a,u], JuMP.getvalue(tp0.assign[a,u,1]))
    end
end

# fixes population a to be un-admixed
function unmix(tp::Union{NodeTreeProblem,TreeProblem}, a::Int)
    JuMP.@constraint(tp.model, 
        [n=getleaves(tp.bt),l=2:tp.nlevels],
        tp.assign[a,n,l] == tp.assign[a,n,l-1])
end

# fixes pop to be un-admixed
function unmix(tp::Union{NodeTreeProblem,TreeProblem}, pop::String)
    a = findfirst(pd.pops .== pop)
    if a == nothing
        @error("Invalid population $(pop)")
    else
        unmix(tp, a)
    end
end

# fixes population a to node u
function fix(tp::Union{NodeTreeProblem,TreeProblem}, a::Int, u::Int)
    JuMP.@constraint(tp.model, 
        [l=1:tp.nlevels],
        tp.assign[a,u,l] == 1)
end

# fixes pop to node u
function fix(tp::Union{NodeTreeProblem,TreeProblem}, pop::String, u::Int)
    a = findfirst(pd.pops .== pop)
    if a == nothing
        @error("Invalid population $(pop)")
    else
        fix(tp, a, u)
    end
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
