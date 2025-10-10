"""
```julia
TurbulenceTypes
```

Module for composite types needed for turbulence transport.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)
"""
module TurbulenceTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes
using ...PinCFlow

include("TurbulencePredictands.jl")
include("TurbulenceIncrements.jl")
include("TurbulenceAuxiliaries.jl")
include("TurbulenceReconstructions.jl")
include("TurbulenceFluxes.jl")
include("Turbulence.jl")

export TurbulencePredictands,
    TurbulenceIncrements,
    TurbulenceAuxiliaries,
    TurbulenceReconstructions,
    TurbulenceFluxes,
    Turbulence
end
