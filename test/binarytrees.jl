module BinaryTreeTests

    using PhylogeneticTrees
    @static if VERSION < v"0.7.0-DEV.2005"
        using Base.Test
    else
        using Test
    end
    
    println("=================================")
    println("Testing BinaryTree data structure")
    println("=================================")

    @testset "Tree construction" begin
        bt4 = PhylogeneticTrees.BinaryTree(4)
        @test bt4.depth == 4
        @test bt4.nnodes == 31
        @test length(bt4.edges) == 30
        @test !PhylogeneticTrees.validnode(bt4,0)
        @test PhylogeneticTrees.validnode(bt4,31)
        @test PhylogeneticTrees.validlayer(bt4,0)
        @test PhylogeneticTrees.validlayer(bt4,4)
        @test PhylogeneticTrees.getlayer(bt4,1) == 0
        @test PhylogeneticTrees.getlayer(bt4,31) == 4
        @test PhylogeneticTrees.getnodes(bt4,0) == 1:1
        @test PhylogeneticTrees.getnodes(bt4,4) == 16:31
        @test all(PhylogeneticTrees.getchildren(bt4,1) .== [2,3])
        @test length(PhylogeneticTrees.getdescendants(bt4,1)) == 30
        @test all(PhylogeneticTrees.getdescendants(bt4,15) .== [30,31])
        @test length(PhylogeneticTrees.getdescendants(bt4,31)) == 0
        @test PhylogeneticTrees.getparent(bt4,31) == 15
    end

    @testset "Path edges" begin 
        bt1 = PhylogeneticTrees.BinaryTree(1)
        @test bt1.pathedges[(1,2)] == [(1, 2)]
        @test all(bt1.pathedges[(2,3)] .== [(1, 2), (1,3)])
        @test length(bt1.pathedges[(1,1)]) == 0
    end

end