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
using ...PinCFlow

include("TracerPredictands.jl")
include("TracerIncrements.jl")
include("TracerAuxiliaries.jl")
include("TracerReconstructions.jl")
include("TracerFluxes.jl")
include("TracerWKBImpact.jl")
include("TracerForcings.jl")
include("Tracer.jl")

export TracerPredictands,
    TracerIncrements,
    TracerAuxiliaries,
    TracerReconstructions,
    TracerFluxes,
    TracerGWForcings,
    TracerWKBImpact,
    Tracer
end
