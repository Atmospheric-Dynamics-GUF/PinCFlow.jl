"""
```julia
Examples
```

Module that defines PinCFlow.jl's example functions.

# Dependencies

  - [`PinCFlow.Examples.WavePacketTools`](@ref)

  - [`PinCFlow.Integration`](@ref)

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.PinCFlow`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.Examples.cold_bubble`](@ref)

  - [`PinCFlow.Examples.hot_bubble`](@ref)

  - [`PinCFlow.Examples.mountain_wave`](@ref)

  - [`PinCFlow.Examples.periodic_hill`](@ref)

  - [`PinCFlow.Examples.vortex`](@ref)

  - [`PinCFlow.Examples.wave_packet`](@ref)

  - [`PinCFlow.Examples.wkb_mountain_wave`](@ref)

  - [`PinCFlow.Examples.wkb_wave_packet`](@ref)
"""
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

export cold_bubble,
    hot_bubble,
    mountain_wave,
    periodic_hill,
    vortex,
    wave_packet,
    wkb_mountain_wave,
    wkb_wave_packet

end
