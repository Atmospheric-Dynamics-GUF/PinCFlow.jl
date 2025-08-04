"""
```julia
IceTypes
```
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
