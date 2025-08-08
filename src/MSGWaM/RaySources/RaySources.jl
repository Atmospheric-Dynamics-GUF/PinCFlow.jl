"""
```julia
RaySources
```

Module for ray-volume sources.

# See also

- [`PinCFlow.Types`](@ref)
- [`PinCFlow.MSGWaM.RayOperations`](@ref)
"""
module RaySources

using ...Types
using ..RayOperations

include("activate_orographic_source!.jl")
include("compute_orographic_mode.jl")

export activate_orographic_source!

end
