include("TestTools/TestTools.jl")

using Test
using .TestTools

submit_directory = "../examples/submit/"

@time @testset verbose = true showtiming = true "PinCFlow tests" begin
    include("test_mountain_wave.jl")
    include("test_periodic_hill.jl")
    include("test_wkb_mountain_wave.jl")
end
