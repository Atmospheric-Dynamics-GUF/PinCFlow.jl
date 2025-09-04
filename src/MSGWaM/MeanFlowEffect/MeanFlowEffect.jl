"""
```julia
MeanFlowEffect
```

Module for computing the mean-flow effect of gravity waves.

Provides functions that compute mean-flow tendencies by integrating ray-volume properties in spectral space and mapping the result to physical grid cells. Also provides two filters for smoothing the tendencies, as well as a simple blocked-layer scheme that includes a blocked-flow drag in mountain-wave simulations.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate`](@ref)
"""
module MeanFlowEffect

using LinearAlgebra
using ...Types
using ...Boundaries
using ..Interpolation
using ..RayUpdate

"""
```julia 
UCHI <: AbstractVariable 
```

Singleton for dispatch to calculation of zonal gravity-wave-tracer fluxes.
"""
struct UCHI <: AbstractVariable end

"""
```julia 
VCHI <: AbstractVariable 
```

Singleton for dispatch to calculation of meridional gravity-wave-tracer fluxes.
"""
struct VCHI <: AbstractVariable end

"""
```julia 
WCHI <: AbstractVariable 
```

Singleton for dispatch to calculation of vertical gravity-wave-tracer fluxes.
"""
struct WCHI <: AbstractVariable end

include("compute_leading_order_tracer_fluxes!.jl")
include("leading_order_tracer_fluxes.jl")
include("compute_leading_order_tracer_forcing!.jl")
include("set_tracer_fields_zero!.jl")
include("apply_blocked_layer_scheme!.jl")
include("apply_shapiro_filter!.jl")
include("compute_gw_integrals!.jl")
include("compute_gw_tendencies!.jl")
include("compute_horizontal_cell_indices.jl")
include("compute_mean_flow_effect!.jl")
include("smooth_gw_tendencies!.jl")

export compute_mean_flow_effect!

end
