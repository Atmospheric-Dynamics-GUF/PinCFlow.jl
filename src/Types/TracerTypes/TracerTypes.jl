"""
```julia
TracerTypes
```

Module for composite types needed for tracer transport.

# Dependencies

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)

# Exports

  - [`PinCFlow.Types.TracerTypes.Tracer`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerFluxes`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerIncrements`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerPredictands`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerReconstructions`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerWKBIntegrals`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerWKBTendencies`](@ref)
"""
module TracerTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes
using ...Macros

include("TracerPredictands.jl")
include("TracerIncrements.jl")
include("TracerReconstructions.jl")
include("TracerFluxes.jl")
include("TracerWKBIntegrals.jl")
include("TracerWKBTendencies.jl")
include("Tracer.jl")

export TracerPredictands,
    TracerIncrements,
    TracerReconstructions,
    TracerFluxes,
    TracerWKBIntegrals,
    TracerWKBTendencies,
    Tracer
end
