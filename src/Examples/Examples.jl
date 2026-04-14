module Examples

include("WavePacketTools/WavePacketTools.jl")

using MPI
using HDF5
using .WavePacketTools
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
include("wp-3d.jl")

export cold_bubble,
    hot_bubble,
    mountain_wave,
    periodic_hill,
    vortex,
    wave_packet,
    wkb_mountain_wave,
    wkb_wave_packet,
    wp_3d

end
