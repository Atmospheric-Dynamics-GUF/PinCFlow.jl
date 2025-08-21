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

  - [`PinCFlow.Types.IceTypes`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes`](@ref)
"""
module Types

"""
```julia
AbstractVariable
```

Abstract type for prognostic variables.
"""
abstract type AbstractVariable end

"""
```julia
Rho <: AbstractVariable
```

Singleton that represents density fluctuations predicted with the continuity equation.
"""
struct Rho <: AbstractVariable end

"""
```julia
RhoP <: AbstractVariable
```

Singleton that represents density fluctuations predicted with the auxiliary equation.
"""
struct RhoP <: AbstractVariable end

"""
```julia
U <: AbstractVariable
```

Singleton that represents the zonal wind.
"""
struct U <: AbstractVariable end

"""
```julia
V <: AbstractVariable
```

Singleton that represents the meridional wind.
"""
struct V <: AbstractVariable end

"""
```julia
W <: AbstractVariable
```

Singleton that represents the (transformed) vertical wind.
"""
struct W <: AbstractVariable end

"""
```julia
PiP <: AbstractVariable
```

Singleton that represents the Exner-pressure fluctuations.
"""
struct PiP <: AbstractVariable end

"""
```julia
P <: AbstractVariable
```

Singleton that represents the mass-weighted potential temperature.
"""
struct P <: AbstractVariable end

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
include("IceTypes/IceTypes.jl")
include("TurbulenceTypes/TurbulenceTypes.jl")

using .NamelistTypes
using .FoundationalTypes
using .PoissonTypes
using .VariableTypes
using .WKBTypes
using .TracerTypes
using .IceTypes
using .TurbulenceTypes

include("State.jl")

export AbstractBackground,
    AbstractCoriolisMode,
    AbstractLimiter,
    AbstractVariable,
    AbstractModel,
    AbstractTestCase,
    AbstractBoundaries,
    AbstractSponge,
    AbstractMergeMode,
    AbstractWKBMode,
    AbstractWKBTestCase,
    AbstractWKBFilter,
    AbstractTracer,
    AbstractIce,
    AbstractTurbulence

export Rho,
    RhoP,
    U,
    V,
    W,
    PiP,
    P,
    Explicit,
    Implicit,
    UniformBoussinesq,
    StratifiedBoussinesq,
    Isothermal,
    FPlane,
    MCVariant,
    Boussinesq,
    PseudoIncompressible,
    Compressible,
    MountainWave,
    WKBMountainWave,
    WavePacket,
    PeriodicBoundaries,
    SolidWallBoundaries,
    ExponentialSponge,
    COSMOSponge,
    PolynomialSponge,
    SinusoidalSponge,
    ConstantWaveAction,
    ConstantWaveEnergy,
    SteadyState,
    SingleColumn,
    MultiColumn,
    Box,
    Shapiro

export DomainNamelist,
    OutputNamelist,
    SettingNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    WavePacketNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    IceNamelist,
    TurbulenceNamelist,
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
    BicGStab,
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
    WKBIncrements,
    SurfaceIndices,
    WKB,
    Tracer,
    Ice,
    Turbulence,
    State,
    NoTracer,
    LinearTracer,
    TracerPredictands,
    TracerAuxiliaries,
    TracerIncrements,
    TracerReconstructions,
    TracerFluxes,
    IceOn,
    NoIce,
    IcePredictands,
    IceAuxiliaries,
    IceIncrements,
    IceReconstructions,
    IceFluxes,
    TurbulenceOn,
    NoTurbulence,
    TurbulencePredictands,
    TurbulenceAuxiliaries,
    TurbulenceIncrements,
    TurbulenceReconstructions,
    TurbulenceFluxes

end
