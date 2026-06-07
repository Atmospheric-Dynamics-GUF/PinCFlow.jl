"""
```julia
@ivy(x::Expr)
```

Return the expression `x` with `@inbounds` and `@views` applied to it.

The `@inbounds` macro is applied to the bodies of all function definitions in `x`.

# Arguments

  - `x`: Input expression.
"""
macro ivy end

macro ivy(x::Expr)
    return esc(ivy(x))
end
