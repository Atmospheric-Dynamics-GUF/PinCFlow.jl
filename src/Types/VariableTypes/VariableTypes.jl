module VariableTypes

using ..NamelistTypes
using ..FoundationalTypes

include("set_p.jl")
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
