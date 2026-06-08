"""
```julia
@ivy(x::Expr)
```

Return the expression `x` with all slices turned into views and all bounds checks turned off.

# Arguments

  - `x`: Input expression.

# See also

  - [`PinCFlow.ivy`](@ref)
"""
macro ivy end

macro ivy(x::Expr)
    return esc(ivy(x))
end
