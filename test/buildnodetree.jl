module BuildTree

    using PhylogeneticTrees, JuMP, Gurobi
    @static if VERSION < v"0.7.0-DEV.2005"
        using Base.Test
    else
        using Test
    end

    println("=============================")
    println("Testing node-based tree model")
    println("=============================")

    @testset "Solving a tree problem" begin

    	pd = PhylogeneticTrees.PopulationData(
    		"testdata/f3.Europe6.csv",
    		"testdata/f3.Europe6-covariance.csv"
    		)

        bt = PhylogeneticTrees.BinaryTree(3)
        tp = PhylogeneticTrees.NodeTreeProblem(pd, bt,
            binaryencoding = true, # this slows it down significantly but is just to test 
            solver = GurobiSolver(OutputFlag = 0))
        PhylogeneticTrees.breaksymmetries(tp, rules = [:leftfirst, :alphabetize])
        @time solve(tp.model)
        # 10.053197 seconds (2.07 M allocations: 107.326 MiB, 0.48% gc time)

        leaves = PhylogeneticTrees.getnodes(bt, bt.depth)

        # test integrality
        for a in 1:pd.npop, u in leaves 
            @test isapprox(getvalue(tp.assign[a,u,1]), 1) || isapprox(getvalue(tp.assign[a,u,1]), 0)
        end

        xval = [leaves[findfirst(round.(getvalue(tp.assign[a,:,1])) .> 0)] for a in 1:pd.npop] 

        @test isapprox(round(getvalue(tp.f3formula[1,1,xval[1],xval[1],1,1])), 482)
        @test isapprox(round(getvalue(tp.f3formula[1,2,xval[1],xval[2],1,1])), 33)
        @test isapprox(round(getvalue(tp.f3formula[2,2,xval[2],xval[2],1,1])), 242)
        for a in 1:pd.npop, b in a:pd.npop, u in leaves, v in leaves 
            u == tp.outgroupnode && continue 
            v == tp.outgroupnode && continue 
            xval[a] == u && continue 
            xval[b] == v && continue 
            @test isapprox(getvalue(tp.f3formula[a,b,u,v,1,1]), 0)
        end

        @test isapprox(getobjectivevalue(tp.model), 302.06620671623)

        PhylogeneticTrees.removesolution(tp, getvalue(tp.assign))
        @time solve(tp.model)
        # 9.538941 seconds (28.87 k allocations: 2.222 MiB)
        @test isapprox(getobjectivevalue(tp.model), 302.06620671623) # symmetric solution
        PhylogeneticTrees.removesolution(tp, getvalue(tp.assign))
        @time solve(tp.model)
        @test isapprox(getobjectivevalue(tp.model), 545.6614084929685) # new solution
        # 10.924954 seconds (13.26 k allocations: 1.406 MiB)

    end

end
