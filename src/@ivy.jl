"""
```julia
@ivy(x::Expr)
```

Return the expression `x` with `@inbounds` and `@views` in front of it.
"""
macro ivy end

macro ivy(x::Expr)
    return esc(quote
        @inbounds @views $x
    end)
end
