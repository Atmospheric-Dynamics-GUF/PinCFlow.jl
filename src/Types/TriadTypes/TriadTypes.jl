module TriadTypes

using ..NamelistTypes
using ..FoundationalTypes
using ..VariableTypes
using ...PinCFlow


include("SpectralGrid.jl")
include("KinematicBox.jl")
include("InterpCoef.jl")
include("ResManifold.jl")
include("log_range.jl")
include("compute_edges.jl")
include("TriadTendencies.jl")


export  TriadTendencies, KinematicBox, SpectralGrid

end
