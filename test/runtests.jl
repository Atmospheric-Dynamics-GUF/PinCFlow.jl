using Test
using PinCFlow_dev
using Trixi: CompressibleEulerEquations1D

@testset "Test SemiDiscretization" begin
    grid = LinRange(0.0, 1.0, 10)
    equations = CompressibleEulerEquations1D(1.4)
    surface_flux = nothing
    initial_condition = nothing

    semi_discretization = SemiDiscretization(grid, equations, surface_flux,
                                             initial_condition)

    @test semi_discretization.grid == grid
end

include("test_poisson.jl")
include("test_elixirs.jl")
