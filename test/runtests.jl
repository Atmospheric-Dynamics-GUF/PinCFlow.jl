include("TestTools/TestTools.jl")

using Test
using PinCFlow
using .TestTools

submit_directory = "../examples/submit/"
update_references = false

@testset verbose = true "PinCFlow tests" begin
    include("mountain_wave.jl")
    include("periodic_hill.jl")
    include("vortex.jl")
    include("wkb_mountain_wave.jl")
end
