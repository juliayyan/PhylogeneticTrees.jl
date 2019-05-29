mutable struct PopulationData
	pops::Vector{String}
	npop::Int
	outgroup::Int
	f3::Matrix{Float64}
	cov::Array{Float64,4}
end 

function PopulationData(
	f3file::String,
	covfile::String)

	f3dat  = CSV.read(f3file)
	covdat = CSV.read(covfile)
	
    # save populations
    pops = sort(vcat(unique(f3dat[:Outgroup]),unique(f3dat[:A])))
	npop = length(pops)
    popinds = Dict(zip(pops,1:npop))
    outgroup = findfirst(pops .== f3dat[:Outgroup][1])

    # get f3(Outgroup; A, B) stats
    f3 = zeros(npop, npop)
    for i in 1:DataFrames.nrow(f3dat)
        a = popinds[f3dat[i,:A]]
        b = popinds[f3dat[i,:B]] 
        f3[a,b] = f3dat[i,:f3]
        f3[b,a] = f3dat[i,:f3]
    end

    # get covariance stats
    cov = zeros(npop, npop, npop, npop);
    for i in 1:DataFrames.nrow(covdat)
        a1 = min(popinds[covdat[i,:A1]],popinds[covdat[i,:B1]])
        b1 = max(popinds[covdat[i,:A1]],popinds[covdat[i,:B1]])
        a2 = min(popinds[covdat[i,:A2]],popinds[covdat[i,:B2]])
        b2 = max(popinds[covdat[i,:A2]],popinds[covdat[i,:B2]])
        cov[a1,b1,a2,b2] = covdat[i,:covariance]
        cov[b1,a1,b2,a2] = covdat[i,:covariance]
        cov[a2,b2,a1,b1] = covdat[i,:covariance]
        cov[b2,a2,b1,a1] = covdat[i,:covariance]
    end

    PopulationData(pops,npop,outgroup,f3,cov)
end

function Base.show(io::IO, pd::PopulationData; offset::String="")
    println(io, offset, string(typeof(pd)))
    println(io, offset, "    * npop: $(pd.npop)  (pops: $(pd.pops[1]), ..., $(pd.pops[end]))")
end

