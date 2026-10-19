module Research

include("WavePacketTools/WavePacketTools.jl")

using MPI
using HDF5
using .WavePacketTools
using ..Macros
using ..Types
using ..Integration
using ..PinCFlow

include("boussinesq_periodic.jl")
include("wkb_boussinesq_periodic.jl")

export boussinesq_periodic, wkb_boussinesq_periodic

end
