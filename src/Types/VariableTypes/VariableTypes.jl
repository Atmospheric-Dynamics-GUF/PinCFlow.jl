"""
```julia
VariableTypes
```

Module for composite types needed for the integration in time.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)

  - [`PinCFlow.Types.FoundationalTypes`](@ref)
"""
module VariableTypes

using ..NamelistTypes
using ..FoundationalTypes
using Random
using ...PinCFlow

include("set_p.jl")

include("../../MSGWaM/RaySources/construct_random_wavepackets!.jl")
include("Predictands.jl")
include("Increments.jl")
include("Backups.jl")
include("Auxiliaries.jl")
include("Reconstructions.jl")
include("Fluxes.jl")
include("Variables.jl")

export Predictands,
    Increments, Backups, Auxiliaries, Reconstructions, Fluxes, Variables

end
