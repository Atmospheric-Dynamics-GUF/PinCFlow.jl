module GWTracerOperations

using LinearAlgebra
using ..BlockedLayer
using ..Interpolation
using ..RayOperations
using ...Types
using ...Boundaries
using ...Update
using ...PinCFlow
using ..Smoothing

"""
```julia
UChi
```

Singleton for dispatch to calculation of zonal gravity-wave-tracer fluxes.
"""
struct UChi end

"""
```julia
VChi
```

Singleton for dispatch to calculation of meridional gravity-wave-tracer fluxes.
"""
struct VChi end

"""
```julia
WChi
```

Singleton for dispatch to calculation of vertical gravity-wave-tracer fluxes.
"""
struct WChi end

include("compute_gw_tracer_integrals!.jl")
include("compute_gw_tracer_tendencies!.jl")
include("compute_next_order_tracer_fluxes!.jl")
include("compute_turbulent_tracer_fluxes!.jl")
include("leading_order_tracer_fluxes.jl")
include("set_tracer_fields_zero!.jl")

export compute_gw_tracer_integrals!,
    compute_gw_tracer_tendencies!,
    compute_next_order_tracer_fluxes!,
    compute_turbulent_tracer_fluxes!,
    leading_order_tracer_fluxes,
    set_tracer_fields_zero!
end