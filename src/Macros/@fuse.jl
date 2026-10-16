"""
```julia
@fuse(x::Any)
```

If `x` is a `for` loop, flatten it and apply `@simd ivdep` to it, otherwise apply FastBroadcast.jl's `@..` to it.

# Argument

  - `x`: Input object.
"""
macro fuse end

macro fuse(x::Any)
    if x isa Expr && x.head === :for
        (loop, body) = x.args

        if loop.head === :block && all(arg.head === :(=) for arg in loop.args)
            index_names = Tuple(arg.args[1] for arg in loop.args[end:(-1):1])
            ranges = Tuple(arg.args[2] for arg in loop.args[end:(-1):1])

            @gensym index indices

            return esc(
                quote
                    let
                        $indices = CartesianIndices(($(ranges...),))

                        @simd ivdep $(
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
                            )
                        )
                    end
                end,
            )
        elseif loop.head === :(=)
            return esc(quote
                @simd ivdep $x
            end)
        else
            error("Unexpected for-loop format!")
        end
    else
        return esc(
            Expr(
                :macrocall,
                GlobalRef(FastBroadcast, Symbol("@..")),
                __source__,
                x,
            ),
        )
    end
end
