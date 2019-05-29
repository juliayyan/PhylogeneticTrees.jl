mutable struct NodeTreeProblem
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

function NodeTreeProblem(
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
    othernodes = leaves[2:end]
    levels = 1:nlevels

    tree = JuMP.Model(solver=solver)
    if binaryencoding
        @assert nlevels == 1
        JuMP.@variable(tree, assign[1:npop,getleaves(bt),levels] >= 0)
        binaryencodingconstraints(pd, bt, tree, assign)
    else 
        JuMP.@variable(tree, assign[1:npop,getleaves(bt),levels], Bin)
    end
    JuMP.@variable(tree, assign2[a=1:npop,b=a:npop,othernodes,othernodes,levels,levels])
    JuMP.@variable(tree, weight[edges] >= 0)
    JuMP.@variable(tree, f3formula[a=1:npop,b=a:npop,u=othernodes,v=othernodes,levels,levels] >= 0)
    JuMP.@variable(tree, f3err[a=1:npop,b=a:npop])
 
    validtreeconstraints(pd, bt, tree, assign, outgroupnode, nlevels, nmixtures)
    logicalconstraints(pd, bt, tree, assign, assign2, outgroupnode, nlevels)
    errorconstraints(pd, bt, tree, assign2, weight, f3formula, f3err, outgroupnode, nlevels)
    
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

    NodeTreeProblem(pd, bt, outgroupnode, nlevels, tree, assign, assign2, weight, f3formula, f3err)
end

function logicalconstraints(
    pd::PopulationData, 
    bt::BinaryTree,
    tree::JuMP.Model, 
    assign::JuMP.JuMPArray{JuMP.Variable},
    assign2::JuMP.JuMPDict{JuMP.Variable},
    outgroupnode::Int,
    nlevels::Int)

    npop = pd.npop
    othernodes = getleaves(bt)[2:end]
    levels = 1:nlevels
    
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

    npop = pd.npop
    othernodes = getleaves(bt)[2:end]
    bigm = maximum(pd.f3)*2
    levels = 1:nlevels

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
