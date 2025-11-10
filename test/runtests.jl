include("TestTools/TestTools.jl")

using Test
using PinCFlow
using .TestTools

submit_directory = "../examples/submit/"
update_references = false

@testset verbose = true "PinCFlow tests" begin
    include("cold_bubble.jl")
    include("hot_bubble.jl")
    include("mountain_wave.jl")
    include("periodic_hill.jl")
    include("vortex.jl")
    include("wave_packet.jl")
    include("wkb_mountain_wave.jl")
end
