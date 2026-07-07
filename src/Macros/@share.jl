macro share end

macro share(x::Vararg{Expr})
    if x[end].head === :for
        block_size = 64
        thread = true
        vector = true

        # Parse the arguments.
        options = Symbol[]
        reductions = Expr[]
        for xi in x[1:(end - 1)]
            if xi.head === :(=)
                for option in options
                    xi.args[1] === option && error("Duplicate option!")
                end

                if xi.args[1] === :block_size && xi.args[2] isa Integer
                    block_size = xi.args[2]
                elseif xi.args[1] === :thread && xi.args[2] isa Bool
                    thread = xi.args[2]
                elseif xi.args[1] === :vector && xi.args[2] isa Bool
                    vector = xi.args[2]
                else
                    error("Unexpected option!")
                end

                push!(options, xi.args[1])
            elseif xi.head === :tuple
                for reduction in reductions
                    xi === reduction && error("Duplicate reduction!")
                end

                if !(xi.args[1] isa Symbol && xi.args[2] isa Symbol)
                    error("Unexpected reduction!")
                end

                push!(reductions, xi)
            else
                error("Unexpected argument format!")
            end
        end

        (loop, body) = x[end].args

        # Extract the index names and ranges of the loop.
        if loop.head === :block && all(arg.head === :(=) for arg in loop.args)
            index_names = Tuple(arg.args[1] for arg in loop.args[end:(-1):1])
            ranges = Tuple(arg.args[2] for arg in loop.args[end:(-1):1])
        elseif loop.head === :(=)
            index_names = (loop.args[1],)
            ranges = (loop.args[2],)
        else
            error("Unexpected for-loop format!")
        end

        @gensym all_indices block blocks index indices loop_size

        vector_loop = Expr(
            :macrocall,
            vector ? Symbol("@simd") : GlobalRef(Macros, Symbol("@identity")),
            __source__,
            :ivdep,
            :(
                for $index in $indices
                    $(Expr(
                        :block,
                        (
                            :($(index_names[i]) = $index[$i]) for
                            i in eachindex(index_names)
                        )...,
                    ))
                    $body
                end
            ),
        )

        thread_loop = Expr(
            :macrocall,
            thread ? GlobalRef(Polyester, Symbol("@batch")) :
            GlobalRef(Macros, Symbol("@identity")),
            __source__,
            length(reductions) > 0 ?
            :(reduction = $(Expr(:tuple, reductions...))) :
            :(stride = false),
            :(
                for $block in $blocks
                    $indices =
                        $all_indices[(($block - 1) * $block_size + 1):min(
                            $block * $block_size,
                            $loop_size,
                        )]
                    $vector_loop
                end
            ),
        )

        output = esc(quote
            $all_indices = CartesianIndices(($(ranges...),))
            $loop_size = length($all_indices)

            $blocks = 1:cld($loop_size, $block_size)
            $thread_loop
        end)

        println(output)

        return output
    else
        if length(x) > 1
            error("Unexpected options!")
        end

        output = esc(
            Expr(
                :macrocall,
                GlobalRef(FastBroadcast, Symbol("@..")),
                __source__,
                :(thread = $(GlobalRef(FastBroadcast, :Threaded))()),
                x[end],
            ),
        )

        println(output)

        return output
    end
end
