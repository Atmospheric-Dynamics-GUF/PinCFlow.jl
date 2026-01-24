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
AbstractTriad
```

Abstract type for triad interactions.
"""
abstract type AbstractTriad end



"""
```julia
AbstractResonance
```

Abstract type for resonance.
"""
abstract type AbstractResonance end




"""
```julia
NeutralStratification <: AbstractBackground
```

Singleton for a Boussinesq atmosphere with neutral stratification.
"""
struct NeutralStratification <: AbstractBackground end

"""
```julia
StableStratification <: AbstractBackground
```

Singleton for a Boussinesq atmosphere with stable stratification.
"""
struct StableStratification <: AbstractBackground end

"""
```julia
RadiatedBoussinesq <: AbstractBackground
```

Singleton for a Boussinesq atmosphere with non uniform stratification genrated by the radiative cooling of the surface.
"""
struct RadiatedBoussinesq <: AbstractBackground end

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

Singleton for switching off MS-GWaM.
"""
struct NoWKB <: AbstractWKBMode end

"""
```julia
SteadyState <: AbstractWKBMode
```

Singleton for the steady-state approximation in MS-GWaM.
"""
struct SteadyState <: AbstractWKBMode end

"""
```julia
SingleColumn <: AbstractWKBMode
```

Singleton for the single-column approximation in MS-GWaM.
"""
struct SingleColumn <: AbstractWKBMode end

"""
```julia
MultiColumn <: AbstractWKBMode
```

Singleton for the multi-column approximation in MS-GWaM.
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
TriadOn <: AbstractTriad
```

Singleton for model configurations with Triad inetractions.
"""
struct TriadOn <: AbstractTriad end

"""
```julia
NoTriad <: AbstractTriad
```

Singleton for model configurations without Triad inetractions.
"""
struct NoTriad <: AbstractTriad end

"""
```julia
Triad2D <: AbstractTriad
```

Singleton for model configurations with Triad inetractions in 2D setup.
"""
struct Triad2D <: AbstractTriad end

"""
```julia
Triad3DIso <: AbstractTriad
```

Singleton for model configurations with Triad inetractions in 3D, horizontally isotropic and vertical axis symmetric setup.
"""
struct Triad3DIso <: AbstractTriad end

"""
```julia
Sum <: AbstractResonance
```

Singleton for model configurations with the sum resonance interactions.
"""
struct Sum <: AbstractResonance end

"""
```julia
Difference <: AbstractResonance
```

Singleton for model configurations with the difference resonance interactions.
"""
struct Difference <: AbstractResonance end


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
include("TriadNamelist.jl")
include("Namelists.jl")


export AbstractBackground,
    AbstractLimiter,
    AbstractModel,
    AbstractMergeMode,
    AbstractWKBMode,
    AbstractWKBFilter,
    AbstractTracer, 
    AbstractTriad,
    AbstractResonance

export NeutralStratification,
    StableStratification,
    RadiatedBoussinesq,
    Isothermal,
    Isentropic,
    Realistic,
    LapseRates,
    MCVariant,
    Boussinesq,
    PseudoIncompressible,
    Compressible,
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
    TriadOn,
    NoTriad,
    Triad2D,
    Triad3DIso, 
    Sum,
    Difference

export DomainNamelist,
    OutputNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TriadNamelist,
    TracerNamelist,
    Namelists

end
