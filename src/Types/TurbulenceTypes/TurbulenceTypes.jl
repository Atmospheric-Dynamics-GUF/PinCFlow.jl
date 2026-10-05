"""
```julia
TurbulenceTypes
```

Module for composite types needed for turbulence parameterization.

# See also

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)
"""
module TurbulenceTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes
using ...Macros

include("TurbulencePredictands.jl")
include("TurbulenceIncrements.jl")
include("TurbulenceAuxiliaries.jl")
include("TurbulenceReconstructions.jl")
include("TurbulenceFluxes.jl")
include("TurbulenceConstants.jl")
include("TurbulenceWKBIntegrals.jl")
include("TurbulenceWKBTendencies.jl")
include("Turbulence.jl")

export TurbulencePredictands,
    TurbulenceIncrements,
    TurbulenceAuxiliaries,
    TurbulenceReconstructions,
    TurbulenceFluxes,
    TurbulenceConstants,
    TurbulenceWKBIntegrals,
    TurbulenceWKBTendencies,
    Turbulence
end
