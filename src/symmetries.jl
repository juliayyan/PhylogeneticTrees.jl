function breaksymmetries(tp::TreeProblem;
    rules::Vector{Symbol} = [:leftfirst])
    bt = tp.bt
    pd = tp.pd
    function addleftfirstrule(cb)
        xval = round.(JuMP.getvalue(tp.assign))
        for layer in getlayer(bt, tp.outgroupnode):(bt.depth-1),
            u in getnodes(bt, layer)
            children = getchildren(bt, u)
            left = children[1]; right = children[2]
            if sum(xval[:,left]) < sum(xval[:,right])
                JuMP.@lazyconstraint(cb, 
                    sum(tp.assign[:,right]) <= sum(tp.assign[:,left]))
            end
        end
    end
    function addalphabetizerule(cb)
        xval = round.(JuMP.getvalue(tp.assign))
        for layer in getlayer(bt, tp.outgroupnode):(bt.depth-1),
            u in getnodes(bt, layer)
            children = getchildren(bt, u)
            left = children[1]; right = children[2]
            if sum(xval[:,left]) + sum(xval[:,right]) == 2
                a = findfirst(xval[:,left])
                b = findfirst(xval[:,right])
                a <= b && continue
                JuMP.@lazyconstraint(cb, 
                    sum(tp.assign[1:b,right]) + sum(tp.assign[(b+1):pd.npop,left]) <= 1)
            end
        end    
    end
    if in(:leftfirst, rules)
        JuMP.addlazycallback(tp.model, addleftfirstrule)
    end
    if in(:alphabetize, rules)
        JuMP.addlazycallback(tp.model, addalphabetizerule)
    end
end
