module BuildTreeCoded

    using PhylogeneticTrees, JuMP, Gurobi
    @static if VERSION < v"0.7.0-DEV.2005"
        using Base.Test
    else
        using Test
    end

    println("========================================")
    println("Testing tree-building with default model")
    println("========================================")

    @testset "Solving a tree problem" begin

    	pd = PhylogeneticTrees.PopulationData(
    		"testdata/f3.Europe6.csv",
    		"testdata/f3.Europe6-covariance.csv"
    		)

        bt = PhylogeneticTrees.BinaryTree(3)
        tp = PhylogeneticTrees.TreeProblem(pd, bt,
            solver = GurobiSolver(OutputFlag = 0),
            binaryencoding = true)
        @time solve(tp.model)
        # 5.566190 seconds (10.90 M allocations: 532.489 MiB, 5.85% gc time)

        leaves = PhylogeneticTrees.getnodes(bt, bt.depth)

        # test integrality
        for a in 1:pd.npop, u in leaves 
            @test isapprox(getvalue(tp.assign[a,u,1]), 1) || isapprox(getvalue(tp.assign[a,u,1]), 0)
        end
        for (a,b,(u,v)) in keys(tp.countedge)
            val = round(getvalue(tp.countedge[a,b,(u,v),1,1]),digits=3)
            @test isapprox(val,1) || 
                  isapprox(val,0)
        end

        xval = [leaves[findfirst(round.(getvalue(tp.assign[a,:,1])) .> 0)] for a in 1:pd.npop] 
        @test isapprox(round(getvalue(tp.f3formula[1,1])), 478)
        @test isapprox(round(getvalue(tp.f3formula[1,2])), 33)
        @test isapprox(round(getvalue(tp.f3formula[2,2])), 247)
        for a in 1:pd.npop, b in a:pd.npop, u in leaves, v in leaves 
            if round(getvalue(tp.assign[a,u,1]) + getvalue(tp.assign[b,v,1])) == 2
                for (u1,v1) in bt.edges 
                    val = round(getvalue(tp.countedge[a,b,(u1,v1),1,1]), digits=3)
                    if in((u1,v1), intersect(bt.pathedges[u,tp.outgroupnode], bt.pathedges[v,tp.outgroupnode]))
                        @test isapprox(val,1)
                    else 
                        @test isapprox(val,0)
                    end
                end
            end
        end

        @test isapprox(getobjectivevalue(tp.model), 267.7057665969791)

    end

end
