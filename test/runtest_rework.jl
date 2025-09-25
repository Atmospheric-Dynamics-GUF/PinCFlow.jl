using Test

const PinCFlow_TEST = get(ENV, "PinCFlow_TEST","all")

@time @testset verbose=true showtiming=true "PinCFlow.jl tests" begin

	@time if PinCFlow_TEST == "all" ## || TODO: 
		include("test_mountain.jl")	
	end

end #end tests
