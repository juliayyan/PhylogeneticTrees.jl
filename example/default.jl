using PhylogeneticTrees, JuMP, Gurobi
using CSV, DataFrames

paramsdf = try CSV.read(ARGS[1]) catch; CSV.read("params-miqo.csv") end
values = [try parse(Int, paramsdf[i,:value]) catch; paramsdf[i,:value] end for i in 1:nrow(paramsdf)]
params = Dict(zip(paramsdf[:param], values))

env = Gurobi.Env()
admix = (params["granularity"] > 1) && (params["admixture_events"] > 0)

# load data
pd = PhylogeneticTrees.PopulationData(params["mean_file"], params["cov_file"])

bt = PhylogeneticTrees.BinaryTree(params["depth"])
tp = PhylogeneticTrees.TreeProblem(pd, bt,
    solver = GurobiSolver(
        env, 
        LogFile = params["log_file"],
        TimeLimit = params["time_limit"]),
    nlevels = params["granularity"],
    nmixtures = pd.npop + params["admixture_events"])
if admix
    if Bool(params["warm_start"])
        PhylogeneticTrees.warmstartunmixed(tp, 
            solver = GurobiSolver(
                env, 
                LogFile = params["log_file"],
                TimeLimit = params["warm_start_time_limit"]))
    end
    if !ismissing(params["unmixed_pops"])
        unmixpops = split(params["unmixed_pops"])
        for pop in unmixpops
            PhylogeneticTrees.unmix(tp, String(pop))        
        end
    end
end
solve(tp.model)
PhylogeneticTrees.printnodes(tp, outfile = params["output_file"])
