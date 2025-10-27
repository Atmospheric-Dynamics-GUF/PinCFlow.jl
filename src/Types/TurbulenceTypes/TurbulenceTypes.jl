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
using ...PinCFlow

include("TurbulencePredictands.jl")
include("TurbulenceIncrements.jl")
include("TurbulenceAuxiliaries.jl")
include("TurbulenceReconstructions.jl")
include("TurbulenceFluxes.jl")
include("TurbulenceConstants.jl")
include("TurbulenceDiffusionCoefficients.jl")
include("Turbulence.jl")

export TurbulencePredictands,
    TurbulenceIncrements,
    TurbulenceAuxiliaries,
    TurbulenceReconstructions,
    TurbulenceFluxes,
    TurbulenceConstants,
    TurbulenceDiffusionCoefficients,
    Turbulence
end
