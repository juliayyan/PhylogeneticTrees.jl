function printnodes(tp::Union{NodeTreeProblem,TreeProblem})
    for u in getnodes(tp.bt, tp.bt.depth), a in 1:tp.pd.npop
        level = round(sum(JuMP.getvalue(tp.assign[a,u,:])))/tp.nlevels
        level > 0 && println(tp.pd.pops[a], "\t", u, "\t", level)
    end
end

# warning: lots of magic constants here
function printtree(tp::Union{NodeTreeProblem,TreeProblem})

    pd = tp.pd
    bt = tp.bt
    depth = bt.depth

    spaces = Dict(zip(0:depth, [3 for i in 0:depth]))
    for d=(depth-1):-1:0
        spaces[d] = 2*spaces[d+1]+1
    end
    for d=0:depth
        nodes = getnodes(bt,d)
        for n in nodes 
            print(" "^spaces[d])
            if d == depth
                nodeassign = round.(JuMP.getvalue(tp.assign[:,n,1]))
            else 
                nodeassign = 0
            end
            if sum(nodeassign) < 1 
                str = "()"
            else 
                str = pd.pops[findfirst(nodeassign .> 0)][1:2]
            end
            Printf.@printf("%2s", str)
            print(" "^spaces[d])
        end
        println()
        if d < depth 
            for n in nodes 
                print(" "^(spaces[d]-2), " /"," "^2,"\\ "," "^(spaces[d]-2))
            end
            println()
            for n in nodes 
                weight1 = JuMP.getvalue(tp.weight[(n,getchildren(bt,n)[1])])
                weight2 = JuMP.getvalue(tp.weight[(n,getchildren(bt,n)[2])])
                print(" "^(spaces[d]-3))
                if round(weight1) > 0
                    Printf.@printf("%3.0f", weight1)
                else 
                    Printf.@printf("%3s", " ")
                end
                print(" "^2)
                if round(weight2) > 0
                    Printf.@printf("%3.0f", weight2)
                else 
                    Printf.@printf("%3s", " ")
                end
                print(" "^(spaces[d]-3))
            end
            println()
            for n in nodes 
                print(" "^(spaces[d]-2),"/ "," "^2," \\"," "^(spaces[d]-2))
            end
            println()
        end  
    end

end

function Base.show(io::IO, tp::Union{NodeTreeProblem,TreeProblem}; offset::String="")
    println(io, offset, string(typeof(tp)))
end
