"""
```julia
@ivy(x::Any)
```

Return `x` with all slices turned into views and all bounds checks turned off.

# Arguments

  - `x`: Input object.

# See also

  - [`PinCFlow.Macros.ivy`](@ref)
"""
macro ivy end

macro ivy(x::Any)
    return esc(ivy(x))
end
