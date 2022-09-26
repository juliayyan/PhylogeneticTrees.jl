mutable struct TreeProblem
    pd::PopulationData
    bt::BinaryTree
    outgroupnode::Int
    nlevels::Int
    model::JuMP.Model
    assign::JuMP.JuMPArray{JuMP.Variable}
    countedge::JuMP.Containers.SparseAxisArray{JuMP.Variable}
    weight::JuMP.JuMPArray{JuMP.Variable}
    f3formula::JuMP.Containers.SparseAxisArray{JuMP.GenericAffExpr{Float64,JuMP.Variable}}
    f3err::JuMP.Containers.SparseAxisArray{JuMP.Variable}
end 

function TreeProblem(
    pd::PopulationData, 
    bt::BinaryTree;
    binaryencoding::Bool = false,
    nlevels::Int = 1,
    nmixtures::Int = pd.npop,
    solver = Gurobi.GurobiSolver())
    
    npop = pd.npop
    edges = bt.edges
    leaves = getleaves(bt)
    outgroupnode = leaves[1]
    levels = 1:nlevels

    tree = JuMP.Model(solver=solver)
    if binaryencoding 
        @assert nlevels == 1
        JuMP.@variable(tree, assign[1:npop,getleaves(bt),levels] >= 0)
        binaryencodingconstraints(pd, bt, tree, assign)
    else 
        JuMP.@variable(tree, assign[1:npop,getleaves(bt),levels], Bin)
    end
    JuMP.@variable(tree, weight[edges] >= 0)
    JuMP.@variable(tree, weightaux[a=1:npop,b=a:npop,edges,levels,levels] >= 0)
    JuMP.@variable(tree, countedge[a=1:npop,b=a:npop,edges,levels,levels] >= 0)
    JuMP.@expression(tree, 
        f3formula[a=1:npop,b=a:npop],
        sum(weightaux[a,b,edg,l,m] for edg in edges, l in levels, m in levels)/nlevels^2)
    JuMP.@variable(tree, f3err[a=1:npop,b=a:npop])

    validtreeconstraints(pd, bt, tree, assign, outgroupnode, nlevels, nmixtures)
    countedgeconstraints(pd, bt, tree, assign, countedge, outgroupnode, nlevels)
    errorconstraints(pd, bt, tree, weight, weightaux, countedge, f3formula, f3err, nlevels)

    keytoind = Dict()
    ind = 0
    for key in keys(f3err)
        key[1] == pd.outgroup && continue
        key[2] == pd.outgroup && continue
        ind += 1
        keytoind[key] = ind
    end
    mat = zeros(ind,ind)
    for (a1,b1) in keys(keytoind), (a2,b2) in keys(keytoind)
        mat[keytoind[a1,b1],keytoind[a2,b2]] = pd.cov[a1,b1,a2,b2]
    end
    matinv = inv(mat)
    
    JuMP.@objective(tree, Min, 
        sum(matinv[keytoind[a1,b1],keytoind[a2,b2]]*f3err[a1,b1]*f3err[a2,b2]
            for (a1,b1) in keys(keytoind), (a2,b2) in keys(keytoind)))

    TreeProblem(
        pd, bt, outgroupnode, nlevels,
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
    countedge::JuMP.Containers.SparseAxisArray{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int)
    
    npop = pd.npop
    levels = 1:nlevels

    for (u,v) in bt.edges 
        if in((u,v), bt.pathedges[outgroupnode,1])
            JuMP.@constraint(tree, 
                [a=1:npop, b=a:npop, l=levels, m=levels],
                countedge[a,b,(u,v),l,m] <= 1-sum(assign[a,n,l] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop, l=levels, m=levels],
                countedge[a,b,(u,v),l,m] <= 1-sum(assign[b,n,m] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop, l=levels, m=levels],
                countedge[a,b,(u,v),l,m] >= 1 - sum(assign[a,n,l] + assign[b,n,m] for n in getsubtreeleaves(bt,v)))
        else
            JuMP.@constraint(tree, 
                [a=1:npop, b=a:npop, l=levels, m=levels],
                countedge[a,b,(u,v),l,m] <= sum(assign[a,n,l] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop, l=levels, m=levels],
                countedge[a,b,(u,v),l,m] <= sum(assign[b,n,m] for n in getsubtreeleaves(bt,v)))
            JuMP.@constraint(tree,  
                [a=1:npop, b=a:npop, l=levels, m=levels],
                countedge[a,b,(u,v),l,m] >= sum(assign[a,n,l] + assign[b,n,m] for n in getsubtreeleaves(bt,v)) - 1)
        end
    end

end

function errorconstraints(pd::PopulationData, 
    bt::BinaryTree, 
    tree::JuMP.Model, 
    weight::JuMP.JuMPArray{JuMP.Variable}, 
    weightaux::JuMP.Containers.SparseAxisArray{JuMP.Variable}, 
    countedge::JuMP.Containers.SparseAxisArray{JuMP.Variable}, 
    f3formula::JuMP.Containers.SparseAxisArray{JuMP.GenericAffExpr{Float64,JuMP.Variable}}, 
    f3err::JuMP.Containers.SparseAxisArray{JuMP.Variable},
    nlevels::Int)

    npop = pd.npop 
    bigm = maximum(pd.f3)*2
    levels = 1:nlevels

    # set weight bilinear terms
    JuMP.@constraint(tree, 
        [a=1:npop,b=a:npop,edg=bt.edges, l=levels, m=levels],
        weightaux[a,b,edg,l,m] <= bigm*countedge[a,b,edg,l,m])
    JuMP.@constraint(tree, 
        [a=1:npop,b=a:npop,edg=bt.edges, l=levels, m=levels],
        weightaux[a,b,edg,l,m] <= weight[edg])
    JuMP.@constraint(tree, 
        [a=1:npop,b=a:npop,edg=bt.edges, l=levels, m=levels],
        weightaux[a,b,edg,l,m] >= weight[edg] + bigm*countedge[a,b,edg,l,m] - bigm)

    # set error terms
    JuMP.@constraint(tree, [a=1:npop,b=a:npop],
        f3err[a,b] == pd.f3[a,b] - f3formula[a,b])
    
end
