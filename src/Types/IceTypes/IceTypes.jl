"""
```julia
IceTypes
```

Module for composite types needed by the ice-physics scheme.

# See also

- [`PinCFlow.Types.NamelistTypes`](@ref)
- [`PinCFlow.Types.FoundationalTypes`](@ref)
- [`PinCFlow.Types.VariableTypes`](@ref)
"""
module IceTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes

include("IcePredictands.jl")
include("IceTendencies.jl")
include("IceAuxiliaries.jl")
include("IceReconstructions.jl")
include("IceFluxes.jl")
include("Ice.jl")

export IcePredictands,
    IceTendencies, IceAuxiliaries, IceReconstructions, IceFluxes, Ice
end
