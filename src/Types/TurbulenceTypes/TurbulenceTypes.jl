"""
```julia
TurbulenceTypes
```

Module for composite types needed by the turbulence scheme.

# See also

- [`PinCFlow.Types.NamelistTypes`](@ref)
- [`PinCFlow.Types.FoundationalTypes`](@ref)
- [`PinCFlow.Types.VariableTypes`](@ref)
"""
module TurbulenceTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes

include("TurbulencePredictands.jl")
include("TurbulenceTendencies.jl")
include("TurbulenceAuxiliaries.jl")
include("TurbulenceReconstructions.jl")
include("TurbulenceFluxes.jl")
include("Turbulence.jl")

export TurbulencePredictands,
    TurbulenceTendencies,
    TurbulenceAuxiliaries,
    TurbulenceReconstructions,
    TurbulenceFluxes,
    Turbulence
end
