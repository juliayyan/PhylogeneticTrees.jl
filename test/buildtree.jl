module BuildTree

    using PhylogeneticTrees, JuMP, Gurobi
    @static if VERSION < v"0.7.0-DEV.2005"
        using Base.Test
    else
        using Test
    end

    println("=====================")
    println("Testing tree-building")
    println("=====================")

    @testset "Solving a tree problem" begin

    	pd = PhylogeneticTrees.PopulationData(
    		"testdata/f3.Europe6.csv",
    		"testdata/f3.Europe6-covariance.csv"
    		)

        bt = PhylogeneticTrees.BinaryTree(3)
        tp = PhylogeneticTrees.TreeProblem(pd, bt, 
            solver = GurobiSolver(OutputFlag = 0))
        PhylogeneticTrees.breaksymmetries(tp, rules = [:leftfirst, :alphabetize])
        @time solve(tp.model)
        #   1.312929 seconds (6.95 k allocations: 3.191 MiB)

        leaves = PhylogeneticTrees.getnodes(bt, bt.depth)
        xval = [leaves[findfirst(round.(getvalue(tp.assign[a,:])))] for a in 1:pd.npop] 
        
        @test isapprox(round(getvalue(tp.f3formula[1,1,xval[1],xval[1]])), 482)
        @test isapprox(round(getvalue(tp.f3formula[1,2,xval[1],xval[2]])), 33)
        @test isapprox(round(getvalue(tp.f3formula[2,2,xval[2],xval[2]])), 242)
        for a in 1:pd.npop, b in a:pd.npop, u in leaves, v in leaves 
            u == tp.outgroupnode && continue 
            v == tp.outgroupnode && continue 
            xval[a] == u && continue 
            xval[b] == v && continue 
            @test isapprox(getvalue(tp.f3formula[a,b,u,v]), 0)
        end

        @test isapprox(getobjectivevalue(tp.model), 302.06620671623)

    end

end
