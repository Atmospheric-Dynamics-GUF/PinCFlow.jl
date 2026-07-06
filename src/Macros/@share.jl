macro share end

macro share(x::Vararg{Expr})
    if x[end].head === :for
        block_size = 64
        simd = true

        options = Symbol[]
        reductions = Expr[]
        for xi in x[1:(end - 1)]
            if xi.head === :(=)
                for option in options
                    xi.args[1] === option && error("Duplicate option!")
                end

                if xi.args[1] === :block_size && xi.args[2] isa Integer
                    block_size = xi.args[2]
                elseif xi.args[1] === :simd && xi.args[2] isa Bool
                    simd = xi.args[2]
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

        if loop.head === :block && all(arg.head === :(=) for arg in loop.args)
            index_names = Tuple(arg.args[1] for arg in loop.args[end:(-1):1])
            ranges = Tuple(arg.args[2] for arg in loop.args[end:(-1):1])
        elseif loop.head === :(=)
            index_names = (loop.args[1],)
            ranges = (loop.args[2],)
        else
            error("Unexpected for-loop format!")
        end

        reduction_argument =
            length(reductions) > 0 ?
            :(reduction = $(Expr(:tuple, reductions...))) : :(stride = false)

        @gensym all_indices block blocks index indices loop_size

        if simd
            output = esc(
                quote
                    $all_indices = CartesianIndices(($(ranges...),))
                    $loop_size = length($all_indices)

                    $blocks = 1:cld($loop_size, $block_size)
                    $(Expr(
                        :macrocall,
                        GlobalRef(Polyester, Symbol("@batch")),
                        __source__,
                        reduction_argument,
                        :(
                            for $block in $blocks
                                $indices =
                                    $all_indices[(($block - 1) * $block_size + 1):min(
                                        $block * $block_size,
                                        $loop_size,
                                    )]

                                @simd for $index in $indices
                                    $(Expr(
                                        :block,
                                        (
                                            :($(index_names[i]) =
                                                    $index[$i]) for
                                            i in eachindex(index_names)
                                        )...,
                                    ))

                                    $body
                                end
                            end
                        ),
                    ))
                end,
            )

            println(output)

            return output
        else
            output = esc(
                Expr(
                    :macrocall,
                    GlobalRef(Polyester, Symbol("@batch")),
                    __source__,
                    reduction_argument,
                    x[end],
                ),
            )

            println(output)

            return output
        end
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
