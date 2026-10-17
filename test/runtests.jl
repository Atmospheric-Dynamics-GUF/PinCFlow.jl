include("TestTools/TestTools.jl")

using JuliaFormatter
using PinCFlow
using Test
using .TestTools

const update = "--update" in ARGS
const keywords = (prepare_restart = true, visualize = false)

include("test_cold_bubble.jl")
include("test_hot_bubble.jl")
include("test_mountain_wave.jl")
include("test_periodic_hill.jl")
include("test_vortex.jl")
include("test_wave_packet.jl")
include("test_wkb_mountain_wave.jl")
include("test_wkb_wave_packet.jl")

include("test_formatting.jl")

@testset verbose = true "PinCFlow tests" begin
    @testset verbose = true "Examples" begin
        test_cold_bubble()
        test_hot_bubble()
        test_mountain_wave()
        test_periodic_hill()
        test_vortex()
        test_wave_packet()
        test_wkb_mountain_wave()
        test_wkb_wave_packet()
    end

    @testset verbose = true "Quality" begin
        test_formatting()
    end
end
