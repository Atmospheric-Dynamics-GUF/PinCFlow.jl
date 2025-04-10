include("../src/PinCFlow.jl")

using .PinCFlow
using Test
using HDF5

@testset "PinCFlow tests" begin
    @testset "Mountain-wave tests" begin
        include("mountain_wave_tests.jl")
    end
end
