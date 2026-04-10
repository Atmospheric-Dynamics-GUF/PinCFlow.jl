"""
```julia
NamelistTypes
```

Module that contains all namelist types.

Provides constructors that allow setting only specific parameters and using default values for the rest.
"""
module NamelistTypes

using MPI
using FunctionWrappers
using ...PinCFlow

import FunctionWrappers: FunctionWrapper

include("@dispatch_background.jl")
include("@dispatch_filter_order.jl")
include("@dispatch_filter_type.jl")
include("@dispatch_limiter_type.jl")
include("@dispatch_merge_mode.jl")
include("@dispatch_model.jl")
include("@dispatch_tracer_setup.jl")
include("@dispatch_wkb_mode.jl")
include("@dispatch_turbulence_scheme.jl")
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
include("TurbulenceNamelist.jl")
include("Namelists.jl")

export @dispatch_background,
    @dispatch_filter_order,
    @dispatch_filter_type,
    @dispatch_limiter_type,
    @dispatch_merge_mode,
    @dispatch_model,
    @dispatch_tracer_setup,
    @dispatch_wkb_mode,
    @dispatch_turbulence_scheme,
    @dispatch

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
