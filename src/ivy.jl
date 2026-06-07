"""
```julia
ivy(x::Any; root::Bool = true)::Any
```

If `x` is an expression, climb its AST and apply `@inbounds` and `@views` to the bodies of scope blocks and macro calls that either of them is not able to pass into, apply `@inbounds` and `@views` to `x` itself if `root == true`, and return it.

# Arguments

  - `x`: Input expression

  - `root`: Switch for applying `@inbounds` and `@views` to `x` itself.
"""
function ivy end

function ivy(x::Any; root::Bool = true)::Any
    if isa(x, Expr)
        if x.head === :macrocall && length(x.args) > 2
            for index in 3:length(x.args)
                if isa(x.args[index], Expr)
                    x.args[index] = quote
                        @inbounds @views $(ivy(x.args[index]; root = false))
                    end
                end
            end
        elseif x.head === :module
            x.args[3] = quote
                @inbounds @views $(ivy(x.args[3]; root = false))
            end
        elseif x.head === :macro ||
               x.head === :function ||
               x.head === :-> ||
               (
                   x.head === :(=) &&
                   isa(x.args[1], Expr) &&
                   x.args[1].head === :call
               ) ||
               x.head === :let
            x.args[2] = quote
                @inbounds @views $(ivy(x.args[2]; root = false))
            end
        elseif x.head === :comprehension || x.head === :generator
            x.args[1] = quote
                @inbounds @views $(ivy(x.args[1]; root = false))
            end
        else
            for index in eachindex(x.args)
                x.args[index] = ivy(x.args[index]; root = false)
            end
        end
    end

    if root
        return quote
            @inbounds @views $x
        end
    else
        return x
    end
end
