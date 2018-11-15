function breaksymmetries(tp::TreeProblem;
    rules::Vector{Symbol} = [:leftfirst])
    bt = tp.bt
    pd = tp.pd
    function addleftfirstrule(cb)
        xval = JuMP.getvalue(tp.assign)
        for u in getnodes(bt, bt.depth-1)
            children = getchildren(bt, u)
            left = children[1]; right = children[2]
            if round(sum(xval[:,left,1])) < round(sum(xval[:,right,1]))
                JuMP.@lazyconstraint(cb, 
                    sum(tp.assign[:,right,1]) <= sum(tp.assign[:,left,1]))
            end
        end
    end
    function addalphabetizerule(cb)
        xval = JuMP.getvalue(tp.assign)
        for u in getnodes(bt, bt.depth-1)
            children = getchildren(bt, u)
            in(tp.outgroupnode, children) && continue
            left = children[1]; right = children[2]
            if round(sum(xval[:,left,1])) + round(sum(xval[:,right,1])) == 2
                a = findfirst(xval[:,left,1])
                b = findfirst(xval[:,right,1])
                a <= b && continue
                JuMP.@lazyconstraint(cb, 
                    sum(tp.assign[1:b,right,1]) + sum(tp.assign[(b+1):pd.npop,left,1]) <= 1)
            end
        end    
    end
    if in(:leftfirst, rules)
        tp.nlevels > 1 && warning("Not tested for admixture")
        JuMP.addlazycallback(tp.model, addleftfirstrule)
    end
    if in(:alphabetize, rules)
        tp.nlevels > 1 && warning("Not tested for admixture")
        JuMP.addlazycallback(tp.model, addalphabetizerule)
    end
end

function removesolution(
    tp::TreeProblem,
    solution::JuMP.JuMPArray{Float64})
    expr = 0
    for a in 1:tp.pd.npop, n in getnodes(tp.bt, tp.bt.depth), l in 1:tp.nlevels
        if solution[a,n,l] > 0
            expr += tp.assign[a,n,l]
        else 
            expr += (1-tp.assign[a,n,l])
        end
    end
    JuMP.@constraint(tp.model, expr <= JuMP.getvalue(expr)-1)
end
