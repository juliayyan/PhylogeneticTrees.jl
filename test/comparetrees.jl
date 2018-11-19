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

    @testset "" begin

    	

    end

end
