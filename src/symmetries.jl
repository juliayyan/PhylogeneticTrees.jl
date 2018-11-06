function breaksymmetries(tp::TreeProblem)
    bt = tp.bt
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
    JuMP.addlazycallback(tp.model, addleftfirstrule)
end
