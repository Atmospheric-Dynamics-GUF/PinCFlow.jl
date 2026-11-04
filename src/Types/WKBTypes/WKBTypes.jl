"""
```julia
WKBTypes
```

Module that contains a collection of types for WKB ray tracing calculations including ray data structures, surface indices, integrals, tendencies, and increments.

# Dependencies

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)

# Exports

  - [`PinCFlow.Types.WKBTypes.MergedRays`](@ref)

  - [`PinCFlow.Types.WKBTypes.Rays`](@ref)

  - [`PinCFlow.Types.WKBTypes.SurfaceIndices`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKB`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBIncrements`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBIntegrals`](@ref)

  - [`PinCFlow.Types.WKBTypes.WKBTendencies`](@ref)
"""
module WKBTypes

import Base: getproperty

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes
using ...Macros

include("WKBIntegrals.jl")
include("WKBTendencies.jl")
include("Rays.jl")
include("MergedRays.jl")
include("WKBIncrements.jl")
include("SurfaceIndices.jl")
include("Spectrum.jl")
include("ElasticModeSelection.jl")
include("WKB.jl")

include("getproperty.jl")

export WKBIntegrals,
    WKBTendencies, Rays, MergedRays, WKBIncrements, SurfaceIndices, WKB

end
