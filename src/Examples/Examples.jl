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

export bhat,
    n2,
    omega,
    phi,
    pihat,
    rhobar,
    thetabar,
    uhat,
    vhat,
    wave_action_density,
    what

export cold_bubble,
    hot_bubble,
    mountain_wave,
    periodic_hill,
    vortex,
    wave_packet,
    wkb_mountain_wave,
    wkb_wave_packet

end
