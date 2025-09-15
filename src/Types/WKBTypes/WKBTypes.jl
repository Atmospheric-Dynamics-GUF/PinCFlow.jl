"""
```julia
WKBTypes
```

Module that contains a collection of types for WKB ray tracing calculations including ray data structures, surface indices, integrals, tendencies, and increments.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)
"""
module WKBTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes
using ...PinCFlow

include("WKBIntegrals.jl")
include("WKBTendencies.jl")
include("Rays.jl")
include("MergedRays.jl")
include("WKBIncrements.jl")
include("SurfaceIndices.jl")
include("WKB.jl")

export WKBIntegrals,
    WKBTendencies, Rays, MergedRays, WKBIncrements, SurfaceIndices, WKB

end
