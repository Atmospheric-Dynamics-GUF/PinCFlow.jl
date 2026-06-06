function deep_inbounds end

function deep_inbounds(
    x::Any;
    current_lnn::Union{Nothing, LineNumberNode} = nothing,
)
    if isa(x, Expr)
        # Update the current LineNumberNode.
        for arg in x.args
            if isa(arg, LineNumberNode)
                current_lnn = arg
            end
        end

        if x.head === :function ||
           x.head === :-> ||
           (x.head === :(=) && isa(x.args[1], Expr) && x.args[1].head === :call)
            # If no LineNumberNode has been encountered, use that of this
            # function.
            lnn =
                current_lnn !== nothing ? current_lnn :
                LineNumberNode(@__LINE__, @__FILE__)

            # Insert @inbounds.
            x.args[2] = Expr(
                :macrocall,
                Symbol("@inbounds"),
                lnn,
                deep_inbounds(x.args[2]; current_lnn),
            )

            return x
        else
            for (index, arg) in enumerate(x.args)
                x.args[index] = deep_inbounds(arg; current_lnn)
            end
        end
    end

    return x
end
