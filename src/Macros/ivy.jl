"""
```julia
ivy(x::Any; root::Bool = true)::Any
```

If `x` is an expression, climb its AST and apply `@inbounds` to the bodies of hard-local-scope blocks (i.e., `macro` and `function` definitions, `let` blocks, comprehensions, and generators), apply `@inbounds` and `@views` to `x` itself if `root == true`, and return it.

# Arguments

  - `x`: Input object.

  - `root`: Switch for applying `@inbounds` and `@views` to `x` itself.
"""
function ivy end

function ivy(x::Any; root::Bool = true)::Any
    if x isa Expr
        if x.head === :macro ||
           x.head === :function ||
           x.head === :-> ||
           (
               x.head === :(=) &&
               isa(x.args[1], Expr) &&
               x.args[1].head === :call
           ) ||
           x.head === :let
            x.args[2] = quote
                @inbounds $(ivy(x.args[2]; root = false))
            end
        elseif x.head === :comprehension || x.head === :generator
            x.args[1] = quote
                @inbounds $(ivy(x.args[1]; root = false))
            end
        else
            for index in eachindex(x.args)
                x.args[index] = ivy(x.args[index]; root = false)
            end
        end
    end

    if root
        return quote
            @views @inbounds $x
        end
    else
        return x
    end
end
