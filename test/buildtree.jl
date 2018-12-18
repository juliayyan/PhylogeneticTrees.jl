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
        # 4.973050 seconds (10.33 M allocations: 504.503 MiB, 6.76% gc time)

        leaves = PhylogeneticTrees.getnodes(bt, bt.depth)

        # test integrality
        for a in 1:pd.npop, u in leaves 
            @test isapprox(getvalue(tp.assign[a,u,1]), 1) || isapprox(getvalue(tp.assign[a,u,1]), 0)
        end
        for (a,b,(u,v)) in keys(tp.countedge)
            @test isapprox(getvalue(tp.countedge[a,b,(u,v),1,1]),1) || 
                  isapprox(getvalue(tp.countedge[a,b,(u,v),1,1]),0)
        end

        xval = [leaves[findfirst(round.(getvalue(tp.assign[a,:,1])) .> 0)] for a in 1:pd.npop] 
        @test isapprox(round(getvalue(tp.f3formula[1,1])), 482)
        @test isapprox(round(getvalue(tp.f3formula[1,2])), 33)
        @test isapprox(round(getvalue(tp.f3formula[2,2])), 242)
        for a in 1:pd.npop, b in a:pd.npop, u in leaves, v in leaves 
            if round(getvalue(tp.assign[a,u,1]) + getvalue(tp.assign[b,v,1])) == 2
                for (u1,v1) in bt.edges 
                    if in((u1,v1), intersect(bt.pathedges[u,tp.outgroupnode], bt.pathedges[v,tp.outgroupnode]))
                        @test isapprox(getvalue(tp.countedge[a,b,(u1,v1),1,1]),1)
                    else 
                        @test isapprox(getvalue(tp.countedge[a,b,(u1,v1),1,1]),0)
                    end
                end
            end
        end        

        @test isapprox(getobjectivevalue(tp.model), 302.06620671623)

    end

end
