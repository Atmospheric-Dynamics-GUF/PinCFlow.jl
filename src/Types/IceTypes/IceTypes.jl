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

include("IceConstants.jl")
include("../../Integration/iceroutines.jl")
include("IcePredictands.jl")
include("IceIncrements.jl")
include("IceAuxiliaries.jl")
include("IceReconstructions.jl")
include("IceFluxes.jl")
include("IceSource.jl")
include("SubGrid.jl")
include("GW.jl")
include("SgsGW.jl")
include("SgsPredictands.jl")
include("Ice.jl")

export IcePredictands,
    IceIncrements, IceAuxiliaries, IceReconstructions, IceFluxes,
    IceSource, IceConstants, Ice, GW,
    SgsGW, SgsPredictands, SubGrid

export psat_ice, sat_ratio, dot_qv, dot_n 

end
