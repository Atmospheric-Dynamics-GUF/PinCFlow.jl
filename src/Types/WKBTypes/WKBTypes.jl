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

include("GWIntegrals.jl")
include("GWTendencies.jl")
include("Rays.jl")
include("Increments.jl")
include("SurfaceIndices.jl")
include("WKB.jl")

export GWIntegrals, GWTendencies, Rays, Increments, SurfaceIndices, WKB

end
