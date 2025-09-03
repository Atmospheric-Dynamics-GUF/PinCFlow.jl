"""
```julia
TracerTypes
```

Module for composite types needed for tracer transport.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)
"""
module TracerTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes

include("TracerPredictands.jl")
include("TracerIncrements.jl")
include("TracerAuxiliaries.jl")
include("TracerReconstructions.jl")
include("TracerFluxes.jl")
include("Tracer.jl")

export TracerPredictands,
    TracerIncrements,
    TracerAuxiliaries,
    TracerReconstructions,
    TracerFluxes,
    Tracer
end
