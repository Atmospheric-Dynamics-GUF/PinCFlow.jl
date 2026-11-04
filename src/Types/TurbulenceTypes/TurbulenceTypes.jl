"""
```julia
TurbulenceTypes
```

Module for composite types needed for turbulence parameterization.

# Scopes

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.VariableTypes`](@ref)

# Exports

  - [`PinCFlow.Types.TurbulenceTypes.Turbulence`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceAuxiliaries`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceConstants`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceFluxes`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceIncrements`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulencePredictands`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceReconstructions`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceWKBIntegrals`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceWKBTendencies`](@ref)
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
