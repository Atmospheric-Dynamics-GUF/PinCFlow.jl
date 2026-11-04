"""
```julia
VariableTypes
```

Module for composite types needed for the integration in time.

# Dependencies

  - [`PinCFlow.Types.FoundationalTypes`](@ref)

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

# Exports

  - [`PinCFlow.Types.VariableTypes.Auxiliaries`](@ref)

  - [`PinCFlow.Types.VariableTypes.Backups`](@ref)

  - [`PinCFlow.Types.VariableTypes.Fluxes`](@ref)

  - [`PinCFlow.Types.VariableTypes.Increments`](@ref)

  - [`PinCFlow.Types.VariableTypes.Predictands`](@ref)

  - [`PinCFlow.Types.VariableTypes.Reconstructions`](@ref)

  - [`PinCFlow.Types.VariableTypes.StressTensor`](@ref)

  - [`PinCFlow.Types.VariableTypes.Variables`](@ref)
"""
module VariableTypes

using ..NamelistTypes
using ..FoundationalTypes
using ...Macros

include("set_p.jl")

include("Predictands.jl")
include("Increments.jl")
include("Backups.jl")
include("Auxiliaries.jl")
include("Reconstructions.jl")
include("Fluxes.jl")
include("StressTensor.jl")
include("Variables.jl")

export Predictands,
    Increments,
    Backups,
    Auxiliaries,
    Reconstructions,
    Fluxes,
    StressTensor,
    Variables

end
