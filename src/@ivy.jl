"""
```julia
@ivy(x::Expr)
```

Return the expression `x` with `@inbounds` and `@views` applied to it and all its subexpressions.

# Arguments

  - `x`: Input expression.

# See also

  - [`PinCFlow.ivy`](@ref)
"""
macro ivy end

macro ivy(x::Expr)
    return esc(ivy(x))
end
