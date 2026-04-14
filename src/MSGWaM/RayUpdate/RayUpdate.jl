"""
```julia
RayUpdate
```

Module for the integration of the ray equations.

In addition to ray-volume initialization and propagation, functions for tracking ray volumes on the model grid and controlling their count, as well as a saturation scheme for capturing wave breaking, are provided.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)

  - [`PinCFlow.MSGWaM.RaySources`](@ref)
"""
module RayUpdate

using MPI
using ..BoundaryRays
using ..Interpolation
using ..RayOperations
using ..RaySources
using ..Smoothing
using ...Types
using ...PinCFlow

include("apply_saturation_scheme!.jl")
include("initialize_rays!.jl")
include("merge_rays!.jl")
include("propagate_rays!.jl")
include("shift_rays!.jl")
include("split_rays!.jl")
include("compute_turbulent_damping.jl")
include("compute_q.jl")
include("compute_turbulent_tracer_fluxes!.jl")
include("compute_gw_turbulence_integrals!.jl")

export apply_saturation_scheme!,
    initialize_rays!,
    merge_rays!,
    propagate_rays!,
    shift_rays!,
    split_rays!,
    compute_q,
    compute_turbulent_tracer_fluxes!,
    compute_gw_turbulence_integrals!

end
