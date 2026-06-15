include("TestTools/TestTools.jl")

using Test
using PinCFlow
using .TestTools

const update_references = false

include("test_cold_bubble.jl")
include("test_hot_bubble.jl")
include("test_mountain_wave.jl")
include("test_periodic_hill.jl")
include("test_vortex.jl")
include("test_wave_packet.jl")
include("test_wkb_mountain_wave.jl")
include("test_wkb_wave_packet.jl")

@testset verbose = true "PinCFlow tests" begin
    test_cold_bubble()
    test_hot_bubble()
    test_mountain_wave()
    test_periodic_hill()
    test_vortex()
    test_wave_packet()
    test_wkb_mountain_wave()
    test_wkb_wave_packet()
end
