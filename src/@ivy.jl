"""
```julia
@ivy(x::Expr)
```

Return the expression `x` with `@inbounds` and `@views` applied to it and all its subexpressions.

This macro places `@inbounds` and `@views` in front of expression arguments of any macro call. Care must therefore be taken when combining `@ivy` with other macros.

# Arguments

  - `x`: Input expression.

# See also

  - [`PinCFlow.ivy`](@ref)
"""
macro ivy end

macro ivy(x::Expr)
    return esc(ivy(x))
end
