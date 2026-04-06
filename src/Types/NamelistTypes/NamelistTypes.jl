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

using MPI
using FunctionWrappers
using ...PinCFlow

import FunctionWrappers: FunctionWrapper

include("@dispatch_filter_order.jl")
include("@dispatch_filter_type.jl")
include("@dispatch_limiter_type.jl")
include("@dispatch_merge_mode.jl")
include("@dispatch_tracer_setup.jl")
include("@dispatch_wkb_mode.jl")
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

export @dispatch_filter_order,
    @dispatch_filter_type,
    @dispatch_limiter_type,
    @dispatch_merge_mode,
    @dispatch_tracer_setup,
    @dispatch_wkb_mode,
    @dispatch

export AbstractBackground, AbstractModel

export NeutralStratification,
    StableStratification,
    Isothermal,
    Isentropic,
    Realistic,
    LapseRates,
    Boussinesq,
    PseudoIncompressible,
    Compressible

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
