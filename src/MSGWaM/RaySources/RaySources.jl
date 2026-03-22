"""
```julia
RaySources
```

Module for ray-volume sources.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MSGWaM.BlockedLayer`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)
"""
module RaySources

using ..BlockedLayer
using ..RayOperations
using ...Types
using ...PinCFlow

include("activate_orographic_source!.jl")
include("compute_orographic_modes!.jl")
include("compute_vertical_averages.jl")

export activate_orographic_source!, compute_orographic_modes!

end
