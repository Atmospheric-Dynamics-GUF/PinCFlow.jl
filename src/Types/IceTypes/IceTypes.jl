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
include("IceActiveVars.jl") # module for ice active variables => must be moved to another place!!
include("IceAuxiliaries.jl")
include("IceReconstructions.jl")
include("IceFluxes.jl")
include("IceSource.jl")
include("SubGrid.jl")
#include("GW.jl")
include("SgsGW.jl")
include("SgsPredictands.jl")
include("SgsIncrements.jl")
include("SgsTendencies.jl")
include("SgsAuxiliaries.jl")
include("Ice.jl")

# dictionary to map the IcePredictands to their corresponding IceAuxiliaryíes (needed for writing output)
const PREDICTAND_TO_AUX = Dict(
    :n_hom => :iaux1,
    :q_hom => :iaux2,
    :qv    => :iaux3,
    :n_in  => :iaux4,
    :n_het => :iaux5,
    :q_het => :iaux6
)

export IcePredictands,
    IceIncrements, IceAuxiliaries, IceReconstructions, IceFluxes,
    IceSource, IceConstants, Ice,
    SgsGW, SgsPredictands, SgsIncrements, SgsTendencies, SgsAuxiliaries, SubGrid, IceActiveVars, get_IceActiveVars, ice_active_vars_tuple

export PREDICTAND_TO_AUX

export psat_ice, sat_ratio, dot_qv, dot_n 

end
