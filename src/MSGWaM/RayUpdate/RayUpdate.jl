"""
```julia
RayUpdate
```

Module for the integration of the ray equations.

In addition to ray-volume initialization and propagation, functions for tracking ray volumes on the model grid and controlling their count, as well as a saturation scheme for capturing wave breaking, are provided.

# Scopes

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)

  - [`PinCFlow.MSGWaM.RaySources`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.MSGWaM.RayUpdate.apply_saturation_scheme!`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate.initialize_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate.merge_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate.propagate_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate.shift_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate.split_rays!`](@ref)
"""
module RayUpdate

using MPI
using ..BoundaryRays
using ..Interpolation
using ..RayOperations
using ..RaySources
using ...Macros
using ...Types

include("apply_saturation_scheme!.jl")
include("initialize_rays!.jl")
include("merge_rays!.jl")
include("propagate_rays!.jl")
include("shift_rays!.jl")
include("split_rays!.jl")
include("apply_turbulent_damping!.jl")

export apply_saturation_scheme!,
    initialize_rays!, merge_rays!, propagate_rays!, shift_rays!, split_rays!

end
