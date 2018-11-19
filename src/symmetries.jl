function breaksymmetries(tp::NodeTreeProblem;
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
        JuMP.addlazycallback(tp.model, addleftfirstrule)
    end
    if in(:alphabetize, rules)
        JuMP.addlazycallback(tp.model, addalphabetizerule)
    end
end
