using Test

const pincflow_test = get(ENV, "PinCFlow_TEST", "all")

@time @testset verbose = true showtiming = true "PinCFlow.jl tests" begin
    @time if pincflow_test == "all" ## || TODO: 
        include("test_mountain.jl")
    end
end #end tests
