"""
```julia
NamelistTypes
```

Module that contains all namelist types.

Provides constructors that allow setting only specific parameters and using default values for the rest.
"""
module NamelistTypes

"""
```julia
AbstractBackground
```

Abstract type for atmospheric background configurations.
"""
abstract type AbstractBackground end

"""
```julia
AbstractCoriolisMode
```

Abstract type for approximations of the Coriolis parameter.
"""
abstract type AbstractCoriolisMode end

"""
```julia
AbstractLimiter
```

Abstract type for flux limiters used in the reconstruction.
"""
abstract type AbstractLimiter end

"""
```julia
AbstractModel
```

Abstract type for levels of compressibility.
"""
abstract type AbstractModel end

"""
```julia
AbstractTestCase
```

Abstract type for model test cases.
"""
abstract type AbstractTestCase end

"""
```julia
AbstractBoundaries
```

Abstract type for vertical boundary conditions.
"""
abstract type AbstractBoundaries end

"""
```julia
AbstractSponge
```

Abstract type for sponge configurations.
"""
abstract type AbstractSponge end

"""
```julia
AbstractMergeMode
```

Abstract type for ray-volume merge algorithms.
"""
abstract type AbstractMergeMode end

"""
```julia
AbstractWKBMode
```

Abstract type for approximations in WKB theory.
"""
abstract type AbstractWKBMode end

"""
```julia
AbstractWKBTestCase <: AbstractTestCase
```

Abstract type for WKB test cases.
"""
abstract type AbstractWKBTestCase <: AbstractTestCase end

"""
```julia
AbstractWKBFilter
```

Abstract type for filtering methods applied to mean-flow tendencies.
"""
abstract type AbstractWKBFilter end

"""
```julia
AbstractTracer
```

Abstract type for the inclusion of a tracer.
"""
abstract type AbstractTracer end

"""
```julia
AbstractIce
```

Abstract type for the inclusion of ice physics.
"""
abstract type AbstractIce end

"""
```julia
AbstractTurbulence
```

Abstract type for the inclusion of turbulence physics.
"""
abstract type AbstractTurbulence end

"""
```julia
UniformBoussinesq <: AbstractBackground
```

Singleton for a Boussinesq atmosphere without stratification.
"""
struct UniformBoussinesq <: AbstractBackground end

"""
```julia
StratifiedBoussinesq <: AbstractBackground
```

Singleton for a Boussinesq atmosphere with stratification.
"""
struct StratifiedBoussinesq <: AbstractBackground end

"""
```julia
Isothermal <: AbstractBackground
```

Singleton for an isothermal atmosphere in pseudo-incompressible or compressible mode.
"""
struct Isothermal <: AbstractBackground end

"""
```julia
FPlane <: AbstractCoriolisMode
```

Singleton for the ``f``-plane approximation of the Coriolis parameter.
"""
struct FPlane <: AbstractCoriolisMode end

"""
```julia
MCVariant <: AbstractLimiter
```

Singleton for the MC-Variant limiter function (used in reconstruction).
"""
struct MCVariant <: AbstractLimiter end

"""
```julia
Boussinesq <: AbstractModel
```

Singleton for Boussinesq dynamics.
"""
struct Boussinesq <: AbstractModel end

"""
```julia
PseudoIncompressible <: AbstractModel
```

Singleton for pseudo-incompressible dynamics.
"""
struct PseudoIncompressible <: AbstractModel end

"""
```julia
Compressible <: AbstractModel
```

Singleton for compressible dynamics.
"""
struct Compressible <: AbstractModel end

"""
```julia
MountainWave <: AbstractTestCase
```

Singleton for mountain-wave test cases.
"""
struct MountainWave <: AbstractTestCase end

"""
```julia
WKBMountainWave <: AbstractWKBTestCase
```

Singleton for WKB-mountain-wave test cases.
"""
struct WKBMountainWave <: AbstractWKBTestCase end

"""
```julia
WavePacket <: AbstractTestCase
```

Singleton for wave-packet test cases.
"""
struct WavePacket <: AbstractTestCase end

"""
```julia
PeriodicBoundaries <: AbstractBoundaries
```

Singleton for periodic boundary conditions in the vertical.
"""
struct PeriodicBoundaries <: AbstractBoundaries end

"""
```julia
SolidWallBoundaries <: AbstractBoundaries
```

Singleton for solid-wall boundary conditions in the vertical.
"""
struct SolidWallBoundaries <: AbstractBoundaries end

"""
```julia
ExponentialSponge <: AbstractSponge
```

Singleton for an exponentially increasing Rayleigh damping in the entire Domain.
"""
struct ExponentialSponge <: AbstractSponge end

"""
```julia
COSMOSponge <: AbstractSponge
```

Singleton for a sponge configuration similar to that used in the COSMO model (squared cosine with time-step-dependent maximum).
"""
struct COSMOSponge <: AbstractSponge end

"""
```julia
PolynomialSponge <: AbstractSponge
```

Singleton for a sponge configuration with polynomial profiles.
"""
struct PolynomialSponge <: AbstractSponge end

"""
```julia
SinusoidalSponge <: AbstractSponge
```

Singleton for a sponge configuration with sinusoidal profiles.
"""
struct SinusoidalSponge <: AbstractSponge end

"""
```julia
ConstantWaveAction <: AbstractMergeMode
```

Singleton for the constant-wave-action ray-volume merging algorithm.
"""
struct ConstantWaveAction <: AbstractMergeMode end

"""
```julia
ConstantWaveEnergy <: AbstractMergeMode
```

Singleton for the constant-wave-energy ray-volume merging algorithm.
"""
struct ConstantWaveEnergy <: AbstractMergeMode end

"""
```julia
SteadyState <: AbstractWKBMode
```

Singleton for the steady-state approximation in MSGWaM.
"""
struct SteadyState <: AbstractWKBMode end

"""
```julia
SingleColumn <: AbstractWKBMode
```

Singleton for the single-column approximation in MSGWaM.
"""
struct SingleColumn <: AbstractWKBMode end

"""
```julia
MultiColumn <: AbstractWKBMode
```

Singleton for the multi-column approximation in MSGWaM.
"""
struct MultiColumn <: AbstractWKBMode end

"""
```julia
Box <: AbstractWKBFilter
```

Singleton for a box filter as smoothing method applied to mean-flow tendencies.
"""
struct Box <: AbstractWKBFilter end

"""
```julia
Shapiro <: AbstractWKBFilter
```

Singleton for a Shapiro filter as smoothing method applied to mean-flow tendencies.
"""
struct Shapiro <: AbstractWKBFilter end

"""
```julia
NoTracer <: AbstractTracer
```

Singleton for model configurations without a tracer.
"""
struct NoTracer <: AbstractTracer end

"""
```julia
LinearTracer <: AbstractTracer
```

Singleton for model configurations with an initially linear tracer.
"""
struct LinearTracer <: AbstractTracer end

"""
```julia
NoIce <: AbstractIce
```

Singleton for model configurations without ice physics.
"""
struct NoIce <: AbstractIce end

"""
```julia
IceOn <: AbstractIce
```

Singleton for model configurations with ice physics.
"""
struct IceOn <: AbstractIce end

"""
```julia
NoTurbulence <: AbstractTurbulence
```

Singleton for model configurations without turbulence physics.
"""
struct NoTurbulence <: AbstractTurbulence end

"""
```julia
TurbulenceOn <: AbstractTurbulence
```

Singleton for model configurations with turbulence physics.
"""
struct TurbulenceOn <: AbstractTurbulence end

struct WKBMultipleWavePackets <: AbstractWKBTestCase end

struct MultipleWavePackets <: AbstractTestCase end

abstract type AbstractCloudCover end

struct CloudCoverOn <: AbstractCloudCover end

struct CloudCoverOff <: AbstractCloudCover end

struct CloudCoverOnLargeScaleOn <: AbstractCloudCover end

abstract type AbstractRandomWavePackets end

struct RandomWavePackets <: AbstractRandomWavePackets end

using MPI

include("DomainNamelist.jl")
include("OutputNamelist.jl")
include("SettingNamelist.jl")
include("DiscretizationNamelist.jl")
include("PoissonNamelist.jl")
include("AtmosphereNamelist.jl")
include("WavePacketNamelist.jl")
include("MultiWavePacketNamelist.jl")
include("GridNamelist.jl")
include("SpongeNamelist.jl")
include("WKBNamelist.jl")
include("TracerNamelist.jl")
include("IceNamelist.jl")
include("TurbulenceNamelist.jl")
include("Namelists.jl")

export AbstractBackground,
    AbstractCoriolisMode,
    AbstractLimiter,
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
    AbstractCloudCover,
    AbstractRandomWavePackets,
    AbstractTurbulence

export UniformBoussinesq,
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
    Shapiro,
    NoTracer,
    LinearTracer,
    NoIce,
    IceOn,
    CloudCoverOff,
    CloudCoverOn,  
    CloudCoverOnLargeScaleOn
    RandomWavePackets,
    NoTurbulence,
    TurbulenceOn, 
    MultipleWavePackets,
    WKBMultipleWavePackets

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
    MultiWavePacketNamelist, 
    Namelists

end
