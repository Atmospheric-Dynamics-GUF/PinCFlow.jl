macro share end

macro share(x::Any)
    if x isa Expr && x.head === :for
        (loop, body) = x.args

        if loop.head === :block && all(arg.head === :(=) for arg in loop.args)
            index_names = Tuple(arg.args[1] for arg in loop.args[end:(-1):1])
            ranges = Tuple(arg.args[2] for arg in loop.args[end:(-1):1])
        elseif loop.head === :(=)
            index_names = (loop.args[1],)
            ranges = (loop.args[2],)
        else
            error("Unexpected for-loop format!")
        end

        @gensym index indices

        output = esc(quote
            let
                $indices = CartesianIndices(($(ranges...),))

                @simd ivdep for $index in $indices
                    $(Expr(
                        :block,
                        (
                            :($(index_names[i]) = $index[$i]) for
                            i in eachindex(index_names)
                        )...,
                    ))
                    $body
                end
            end
        end)

        println(output)

        return output
    else
        output = esc(
            Expr(
                :macrocall,
                GlobalRef(FastBroadcast, Symbol("@..")),
                __source__,
                x,
            ),
        )

        println(output)

        return output
    end
end
