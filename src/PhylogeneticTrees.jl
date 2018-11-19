module PhylogeneticTrees

    import JuMP, Gurobi, MathProgBase
    import LightGraphs
    import CSV, Iterators, DataFrames

    include("tree.jl")
    include("popdata.jl")
    include("model-coded.jl")
    include("model.jl")
    include("symmetries.jl")
end
