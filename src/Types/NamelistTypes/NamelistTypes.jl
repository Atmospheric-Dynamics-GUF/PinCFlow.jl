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
BoxFilter <: AbstractWKBFilter
```

Singleton for a box filter as smoothing method applied to mean-flow tendencies.
"""
struct BoxFilter <: AbstractWKBFilter end

"""
```julia
ShapiroFilter <: AbstractWKBFilter
```

Singleton for a Shapiro filter as smoothing method applied to mean-flow tendencies.
"""
struct ShapiroFilter <: AbstractWKBFilter end

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

using MPI
using FunctionWrappers
using ...PinCFlow

import FunctionWrappers: FunctionWrapper

include("@dispatch_limiter_type.jl")
include("@dispatch.jl")

include("DomainNamelist.jl")
include("OutputNamelist.jl")
include("DiscretizationNamelist.jl")
include("PoissonNamelist.jl")
include("AtmosphereNamelist.jl")
include("GridNamelist.jl")
include("SpongeNamelist.jl")
include("WKBNamelist.jl")
include("TracerNamelist.jl")
include("Namelists.jl")

export @dispatch_limiter_type, @dispatch

export AbstractBackground,
    AbstractModel,
    AbstractMergeMode,
    AbstractWKBMode,
    AbstractWKBFilter,
    AbstractTracer

export NeutralStratification,
    StableStratification,
    Isothermal,
    Isentropic,
    Realistic,
    LapseRates,
    Boussinesq,
    PseudoIncompressible,
    Compressible,
    ConstantWaveAction,
    ConstantWaveEnergy,
    NoWKB,
    SteadyState,
    SingleColumn,
    MultiColumn,
    BoxFilter,
    ShapiroFilter,
    NoTracer,
    TracerOn

export DomainNamelist,
    OutputNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    Namelists

end
