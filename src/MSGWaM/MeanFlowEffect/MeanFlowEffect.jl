"""
```julia
MeanFlowEffect
```

Module for computing the mean-flow effect of gravity waves.

Provides functions that compute mean-flow tendencies by integrating ray-volume properties in spectral space and mapping the result to physical grid cells. Also provides two filters for smoothing the tendencies, as well as a simple blocked-layer scheme that includes a blocked-flow drag in mountain-wave simulations.

# See also

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.MSGWaM.BlockedLayer`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate`](@ref)
"""
module MeanFlowEffect

using LinearAlgebra
using ..BlockedLayer
using ..Interpolation
using ..RayUpdate
using ...Macros
using ...Types
using ...Boundaries

"""
```julia 
Tendencies
```

Singleton for dispatch to smoothing of gravity-wave tendencies.
"""
struct Tendencies end

"""
```julia 
Integrals
```

Singleton for dispatch to smoothing of gravity-wave integrals.
"""
struct Integrals end

"""
```julia 
LeadingOrder
```

Singleton for dispatch to the calculation of the first-order gravity-wave integrals.
"""
struct LeadingOrder end

"""
```julia 
NextOrder
```

Singleton for dispatch to the calculation of the second-order gravity-wave integrals.
"""
struct NextOrder end

include("compute_gw_tracer_integrals!.jl")
include("compute_gw_tracer_tendencies!.jl")
include("reset_tracer_fields!.jl")
include("apply_shapiro_filter!.jl")
include("compute_gw_integrals!.jl")
include("compute_gw_tendencies!.jl")
include("compute_horizontal_cell_indices.jl")
include("compute_mean_flow_effect!.jl")
include("smoothing!.jl")

export compute_mean_flow_effect!, compute_gw_integrals!, smoothing!

export Integrals, LeadingOrder, NextOrder

end
