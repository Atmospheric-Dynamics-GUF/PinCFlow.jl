"""
```julia
Types
```

Module for the construction of a single composite type that contains all information on the current model state.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.PoissonTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)

  - [`PinCFlow.Types.WKBTypes`](@ref)

  - [`PinCFlow.Types.TracerTypes`](@ref)
"""
module Types

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

include("NamelistTypes/NamelistTypes.jl")
include("FoundationalTypes/FoundationalTypes.jl")
include("PoissonTypes/PoissonTypes.jl")
include("VariableTypes/VariableTypes.jl")
include("WKBTypes/WKBTypes.jl")
include("TracerTypes/TracerTypes.jl")

using .NamelistTypes
using .FoundationalTypes
using .PoissonTypes
using .VariableTypes
using .WKBTypes
using .TracerTypes
using ..PinCFlow

include("State.jl")

export @dispatch_background,
    @dispatch_filter_order,
    @dispatch_filter_type,
    @dispatch_limiter_type,
    @dispatch_merge_mode,
    @dispatch_model,
    @dispatch_tracer_setup,
    @dispatch_wkb_mode,
    @dispatch

export AbstractPredictand

export Rho, RhoP, U, V, W, PiP, P, Theta, Chi, Explicit, Implicit

export DomainNamelist,
    OutputNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    Namelists,
    Time,
    Constants,
    Domain,
    Grid,
    Atmosphere,
    Sponge,
    Tensor,
    Operator,
    Preconditioner,
    BiCGSTAB,
    Correction,
    Poisson,
    Predictands,
    Increments,
    Backups,
    Auxiliaries,
    Reconstructions,
    Fluxes,
    Variables,
    WKBIntegrals,
    WKBTendencies,
    Rays,
    MergedRays,
    WKBIncrements,
    SurfaceIndices,
    WKB,
    Tracer,
    State,
    TracerPredictands,
    TracerAuxiliaries,
    TracerIncrements,
    TracerReconstructions,
    TracerFluxes,
    TracerWKBIntegrals,
    TracerWKBTendencies

end
