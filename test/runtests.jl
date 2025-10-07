using Test

@time @testset verbose = true showtiming = true "PinCFlow.jl tests" begin
    include("test_periodic_hill.jl")
end #end tests