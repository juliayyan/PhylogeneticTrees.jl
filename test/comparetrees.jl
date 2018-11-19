module CompareTrees

    using PhylogeneticTrees, JuMP, Gurobi
    @static if VERSION < v"0.7.0-DEV.2005"
        using Base.Test
    else
        using Test
    end

    println("=====================")
    println("Comparing both models")
    println("=====================")

    @testset "Comparing admixture models" begin

        pd = PhylogeneticTrees.PopulationData(
            "testdata/f3.Europe6.csv",
            "testdata/f3.Europe6-covariance.csv"
            )

        nlevels = 2
        depth = 4
        bt = PhylogeneticTrees.BinaryTree(depth)

        tp = PhylogeneticTrees.TreeProblem(pd, bt,
            solver = GurobiSolver(OutputFlag = 0, TimeLimit = 1*60),
            nlevels = nlevels)
        PhylogeneticTrees.warmstartunmixed(tp)
        @time solve(tp.model)
        # 60.237210 seconds (53.08 k allocations: 9.018 MiB)

        tp2 = PhylogeneticTrees.NodeTreeProblem(pd, bt,
            nlevels = nlevels, 
            solver = GurobiSolver(OutputFlag = 0))
        for p in 1:pd.npop, u in PhylogeneticTrees.getleaves(bt), l in 1:nlevels 
            JuMP.@constraint(tp2.model, tp2.assign[p,u,l] == getvalue(tp.assign[p,u,l]))
        end
        @time solve(tp2.model)
        # 0.617204 seconds (155 allocations: 47.423 MiB, 24.36% gc time)

        @test isapprox(getobjectivevalue(tp2.model), getobjectivevalue(tp.model))

    end

end
