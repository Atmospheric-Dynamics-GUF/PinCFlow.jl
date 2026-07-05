macro share end

macro share(x::Vararg{Expr})
    if x[end].head === :for
        block_size = 64
        simd = true

        options = Symbol[]
        for xi in x[1:(end - 1)]
            xi.head !== :(=) && error("Unexpected option format!")

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
        end

        (loop, body) = x[end].args

        if loop.head === :block && all(arg.head === :(=) for arg in loop.args)
            index_names = Tuple(arg.args[1] for arg in loop.args[end:(-1):1])
            ranges = Tuple(
                arg.args[2] isa CartesianIndices ? arg.args[2].indices :
                arg.args[2] for arg in loop.args[end:(-1):1]
            )
        elseif loop.head === :(=)
            index_names = (loop.args[1],)
            ranges = (
                loop.args[2] isa CartesianIndices ? loop.args[2].indices :
                loop.args[2],
            )
        else
            error("Unexpected for-loop format!")
        end

        println(index_names)
        println(ranges)

        @gensym all_indices block blocks index indices loop_size

        if simd
            output = quote
                $all_indices = CartesianIndices(($(esc.(ranges)...),))
                $loop_size = length($all_indices)

                $blocks = 1:cld($loop_size, $block_size)
                $(Expr(
                    :macrocall,
                    GlobalRef(Polyester, Symbol("@batch")),
                    __source__,
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
                                        :($(esc(index_names[i])) =
                                                $index[$i]) for i in 1:length(index_names)
                                    )...,
                                ))
                                $(esc(body))
                            end
                        end
                    ),
                ))
            end

            println(output)

            return output
        else
            output = Expr(
                :macrocall,
                GlobalRef(Polyester, Symbol("@batch")),
                __source__,
                esc(x[end]),
            )

            println(output)

            return output
        end
    else
        if length(x) > 1
            error("Unexpected options!")
        end

        output = Expr(
            :macrocall,
            GlobalRef(FastBroadcast, Symbol("@..")),
            __source__,
            :(thread = $(GlobalRef(FastBroadcast, :Threaded))()),
            esc(x[end]),
        )

        println(output)

        return output
    end
end
