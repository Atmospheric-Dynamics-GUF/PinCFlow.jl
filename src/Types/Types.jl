"""
```julia
Types
```

Module for the construction of a single composite type that contains all information on the current model state.

# Dependencies

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.PoissonTypes`](@ref)

  - [`PinCFlow.Types.TracerTypes`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)

  - [`PinCFlow.Types.WKBTypes`](@ref)

# Exports

  - [`PinCFlow.Types.AbstractPredictand`](@ref)

  - [`PinCFlow.Types.Chi`](@ref)

  - [`PinCFlow.Types.Explicit`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Domain`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Grid`](@ref)

  - [`PinCFlow.Types.Implicit`](@ref)

  - [`PinCFlow.Types.NamelistTypes.AtmosphereNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.DiscretizationNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.DomainNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.GridNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.Namelists`](@ref)

  - [`PinCFlow.Types.NamelistTypes.OutputNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.PoissonNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.SpongeNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.TracerNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.TurbulenceNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.WKBNamelist`](@ref)

  - [`PinCFlow.Types.P`](@ref)

  - [`PinCFlow.Types.PiP`](@ref)

  - [`PinCFlow.Types.Rho`](@ref)

  - [`PinCFlow.Types.RhoP`](@ref)

  - [`PinCFlow.Types.State`](@ref)

  - [`PinCFlow.Types.Theta`](@ref)

  - [`PinCFlow.Types.TKE`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerFluxes`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerPredictands`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerReconstructions`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerWKBIntegrals`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerWKBTendencies`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceFluxes`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulencePredictands`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceReconstructions`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceWKBIntegrals`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceWKBTendencies`](@ref)

  - [`PinCFlow.Types.U`](@ref)

  - [`PinCFlow.Types.V`](@ref)

  - [`PinCFlow.Types.VariableTypes.Fluxes`](@ref)

  - [`PinCFlow.Types.VariableTypes.Predictands`](@ref)

  - [`PinCFlow.Types.VariableTypes.Reconstructions`](@ref)

  - [`PinCFlow.Types.VariableTypes.StressTensor`](@ref)

  - [`PinCFlow.Types.W`](@ref)

  - [`PinCFlow.Types.WKBTypes.MergedRays`](@ref)

  - [`PinCFlow.Types.WKBTypes.Rays`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBIncrements`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBIntegrals`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBTendencies`](@ref)

  - [`PinCFlow.Types.X`](@ref)

  - [`PinCFlow.Types.XZ`](@ref)

  - [`PinCFlow.Types.XYZ`](@ref)

  - [`PinCFlow.Types.Y`](@ref)

  - [`PinCFlow.Types.YZ`](@ref)

  - [`PinCFlow.Types.Z`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_background`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_filter_order`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_filter_type`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_limiter_type`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_merge_mode`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_model`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_orographic_flow`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_tracer_setup`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_turbulence_scheme`](@ref)

  - [`PinCFlow.Types.NamelistTypes.@dispatch_wkb_mode`](@ref)
"""
module Types

using MPI

"""
```julia
AbstractPredictand
```

Abstract type for prognostic variables.
"""
abstract type AbstractPredictand end

"""
```julia
Rho <: AbstractPredictand
```

Singleton that represents density fluctuations predicted with the continuity equation.
"""
struct Rho <: AbstractPredictand end

"""
```julia
RhoP <: AbstractPredictand
```

Singleton that represents density fluctuations predicted with the auxiliary equation.
"""
struct RhoP <: AbstractPredictand end

"""
```julia
U <: AbstractPredictand
```

Singleton that represents the zonal wind.
"""
struct U <: AbstractPredictand end

"""
```julia
V <: AbstractPredictand
```

Singleton that represents the meridional wind.
"""
struct V <: AbstractPredictand end

"""
```julia
W <: AbstractPredictand
```

Singleton that represents the (transformed) vertical wind.
"""
struct W <: AbstractPredictand end

"""
```julia
PiP <: AbstractPredictand
```

Singleton that represents the Exner-pressure fluctuations.
"""
struct PiP <: AbstractPredictand end

"""
```julia
P <: AbstractPredictand
```

Singleton that represents the mass-weighted potential temperature.
"""
struct P <: AbstractPredictand end

"""
```julia
Theta
```

Singleton that represents the potential temperature.
"""
struct Theta end

"""
```julia
Chi
```

Singleton that represents the tracer mixing ratio.
"""
struct Chi end

"""
```julia
TKE
```

Singleton that represents the turbulent kinetic energy.
"""
struct TKE end

"""
```julia
Explicit
```

Singleton for explicit integration in time.
"""
struct Explicit end

"""
```julia
Implicit
```

Singleton for implicit integration in time.
"""
struct Implicit end

"""
```julia
X
```

Singleton for dispatch to operations in ``x``-direction.
"""
struct X end

"""
```julia
Y
```

Singleton for dispatch to operations in ``y``-direction.
"""
struct Y end

"""
```julia
Z
```

Singleton for dispatch to operations in ``z``-direction.
"""
struct Z end

"""
```julia
XZ
```

Singleton for dispatch to operations in ``x``- and ``z``-direction.
"""
struct XZ end

"""
```julia
YZ
```

Singleton for dispatch to operations in ``y``- and ``z``-direction.
"""
struct YZ end

"""
```julia
XYZ
```

Singleton for dispatch to operations in all directions.
"""
struct XYZ end

include("NamelistTypes/NamelistTypes.jl")
include("FoundationalTypes/FoundationalTypes.jl")
include("PoissonTypes/PoissonTypes.jl")
include("VariableTypes/VariableTypes.jl")
include("WKBTypes/WKBTypes.jl")
include("TracerTypes/TracerTypes.jl")
include("TurbulenceTypes/TurbulenceTypes.jl")

using .NamelistTypes
using .FoundationalTypes
using .PoissonTypes
using .VariableTypes
using .WKBTypes
using .TracerTypes
using .TurbulenceTypes

include("State.jl")

export @dispatch_background,
    @dispatch_filter_order,
    @dispatch_filter_type,
    @dispatch_limiter_type,
    @dispatch_merge_mode,
    @dispatch_model,
    @dispatch_orographic_flow,
    @dispatch_tracer_setup,
    @dispatch_turbulence_scheme,
    @dispatch_wkb_mode

export AbstractPredictand

export Rho,
    RhoP,
    U,
    V,
    W,
    PiP,
    P,
    Theta,
    Chi,
    TKE,
    Explicit,
    Implicit,
    X,
    Y,
    Z,
    XZ,
    YZ,
    XYZ

export DomainNamelist,
    OutputNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    TurbulenceNamelist,
    Namelists,
    Domain,
    Grid,
    Predictands,
    Reconstructions,
    Fluxes,
    StressTensor,
    WKBIntegrals,
    WKBTendencies,
    Rays,
    MergedRays,
    WKBIncrements,
    State,
    TracerPredictands,
    TracerReconstructions,
    TracerFluxes,
    TracerWKBIntegrals,
    TracerWKBTendencies,
    TurbulencePredictands,
    TurbulenceReconstructions,
    TurbulenceFluxes,
    TurbulenceWKBIntegrals,
    TurbulenceWKBTendencies

end
