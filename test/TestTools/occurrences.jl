"""
```julia
occurrences(f::Function, x::Any)::Integer
```

Return the number of occurences of a pattern defined by `f` in an expression `x`.

The function `f` is expected to take a single argument `x` (which may be an expression) and return a `Bool`.

# Arguments

  - `f`: Function that defines the pattern of interest.

  - `x`: Expression in which to count the occurences of the pattern defined by `f`.
"""
function occurrences end

function occurrences(f::Function, x::Any)::Integer
    count = f(x) ? 1 : 0

    if x isa Expr
        for arg in x.args
            count += occurrences(f, arg)
        end
    end

    return count
end
