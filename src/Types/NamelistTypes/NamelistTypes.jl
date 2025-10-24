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
AbstractTurbulence
```

Abstract type for the inclusion of a turbulence parameterization.
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
Isentropic <: AbstractBackground
```

Singleton for an isentropic atmosphere in pseudo-incompressible or compressible mode.
"""
struct Isentropic <: AbstractBackground end

"""
```julia
Realistic <: AbstractBackground
```

Singleton for a realistic atmosphere in pseudo-incompressible or compressible mode (isentropic troposphere and isothermal stratosphere).
"""
struct Realistic <: AbstractBackground end

"""
```julia
LapseRates <: AbstractBackground
```

Singleton for an atmosphere with different lapse rates in the troposphere and stratosphere in pseudo-incompressible or compressible mode.
"""
struct LapseRates <: AbstractBackground end

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
NoWKB <: AbstractWKBMode
```

Singleton for switching off MSGWaM.
"""
struct NoWKB <: AbstractWKBMode end

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
TracerOn <: AbstractTracer
```

Singleton for model configurations with an initially linear tracer.
"""
struct TracerOn <: AbstractTracer end

"""
```julia 
NoTurbulence <: AbstractTurbulence
```

Singleton for model configurations without turbulence parameterization.
"""
struct NoTurbulence <: AbstractTurbulence end 

"""
```julia 
TKEScheme <: AbstractTurbulence
```

Singleton for model configurations with turbulence parameterization using a TKE-Scheme.
"""
struct TKEScheme <: AbstractTurbulence end 

using MPI
using ...PinCFlow

include("DomainNamelist.jl")
include("OutputNamelist.jl")
include("DiscretizationNamelist.jl")
include("PoissonNamelist.jl")
include("AtmosphereNamelist.jl")
include("GridNamelist.jl")
include("SpongeNamelist.jl")
include("WKBNamelist.jl")
include("TracerNamelist.jl")
include("TurbulenceNamelist.jl")
include("Namelists.jl")

export AbstractBackground,
    AbstractLimiter,
    AbstractModel,
    AbstractSponge,
    AbstractMergeMode,
    AbstractWKBMode,
    AbstractWKBFilter,
    AbstractTracer,
    AbstractTurbulence

export UniformBoussinesq,
    StratifiedBoussinesq,
    Isothermal,
    Isentropic,
    Realistic,
    LapseRates,
    MCVariant,
    Boussinesq,
    PseudoIncompressible,
    Compressible,
    ExponentialSponge,
    COSMOSponge,
    PolynomialSponge,
    SinusoidalSponge,
    ConstantWaveAction,
    ConstantWaveEnergy,
    NoWKB,
    SteadyState,
    SingleColumn,
    MultiColumn,
    Box,
    Shapiro,
    NoTracer,
    TracerOn,
    NoTurbulence,
    TKEScheme

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
    Namelists

end
