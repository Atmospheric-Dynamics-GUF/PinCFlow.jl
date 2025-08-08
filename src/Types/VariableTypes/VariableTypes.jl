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

include("Predictands.jl")
include("Tendencies.jl")
include("Backups.jl")
include("Auxiliaries.jl")
include("Reconstructions.jl")
include("Fluxes.jl")
include("Variables.jl")

export Predictands,
    Tendencies, Backups, Auxiliaries, Reconstructions, Fluxes, Variables

end
