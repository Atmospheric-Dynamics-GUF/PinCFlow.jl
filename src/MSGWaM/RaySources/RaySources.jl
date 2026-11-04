"""
```julia
RaySources
```

Module for ray-volume sources.

# Dependencies

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.MSGWaM.BlockedLayer`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref)

  - [`PinCFlow.MSGWaM.RaySources.compute_orographic_modes!`](@ref)
"""
module RaySources

using ..BlockedLayer
using ..RayOperations
using ...Macros
using ...Types

include("activate_orographic_source!.jl")
include("apply_elastic_mode_selection!.jl")
include("compute_orographic_modes!.jl")

export activate_orographic_source!, compute_orographic_modes!

end
