module Examples

include("WavePacketTools/WavePacketTools.jl")

using MPI
using HDF5
using .WavePacketTools
using ..Macros
using ..Types
using ..Integration
using ..PinCFlow

include("cold_bubble.jl")
include("hot_bubble.jl")
include("mountain_wave.jl")
include("periodic_hill.jl")
include("vortex.jl")
include("wave_packet.jl")
include("wkb_mountain_wave.jl")
include("wkb_wave_packet.jl")
include("wave_Boussinesq.jl")
include("wkb_wave_Boussinesq.jl")

export cold_bubble,
    hot_bubble,
    mountain_wave,
    periodic_hill,
    vortex,
    wave_packet,
    wkb_mountain_wave,
    wkb_wave_packet,
    wave_Boussinesq,
    wkb_wave_Boussinesq

end
