function deep_inbounds end

function deep_inbounds(x::Any)::Any
    if isa(x, Expr)
        if x.head === :function ||
           x.head === :-> ||
           (x.head === :(=) && isa(x.args[1], Expr) && x.args[1].head === :call)
            x.args[2] = Expr(
                :macrocall,
                Symbol("@inbounds"),
                LineNumberNode(@__LINE__, @__FILE__),
                deep_inbounds(x.args[2]),
            )
        else
            for (index, arg) in enumerate(x.args)
                x.args[index] = deep_inbounds(arg)
            end
        end
    end

    return x
end
