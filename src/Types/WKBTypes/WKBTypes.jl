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
include("SpectralGrid.jl")
include("KinematicBox.jl")
include("InterpCoef.jl")
include("TriadTendencies.jl")
include("log_range.jl")
include("compute_edges.jl")
include("WKB.jl")

export WKBIntegrals,
    WKBTendencies, Rays, MergedRays, WKBIncrements, SurfaceIndices, TriadTendencies, WKB, KinematicBox, SpectralGrid

end
