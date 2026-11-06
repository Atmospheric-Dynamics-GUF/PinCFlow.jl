"""
```julia
@ivy(x::Any)
```

Return `x` with all slices turned into views and all bounds checks turned off.

# Arguments

  - `x`: Input object.

# Calls

  - [`PinCFlow.Macros.ivy`](@ref)
"""
macro ivy end

macro ivy(x::Any)
    esc(ivy(x))
end
