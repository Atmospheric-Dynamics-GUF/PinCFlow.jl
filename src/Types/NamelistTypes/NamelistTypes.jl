"""
```julia
NamelistTypes
```

Module that contains all namelist types.

Provides constructors that allow setting only specific parameters and using default values for the rest.

# Scopes

  - [`PinCFlow.Macros`](@ref)

# Exports

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
"""
module NamelistTypes

using FunctionWrappers: FunctionWrapper
using ...Macros

include("@dispatch_background.jl")
include("@dispatch_filter_order.jl")
include("@dispatch_filter_type.jl")
include("@dispatch_limiter_type.jl")
include("@dispatch_merge_mode.jl")
include("@dispatch_model.jl")
include("@dispatch_orographic_flow.jl")
include("@dispatch_tracer_setup.jl")
include("@dispatch_turbulence_scheme.jl")
include("@dispatch_wkb_mode.jl")

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
    @dispatch_orographic_flow,
    @dispatch_tracer_setup,
    @dispatch_turbulence_scheme,
    @dispatch_wkb_mode

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
