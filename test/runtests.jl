using Test
using PinCFlow_dev
using Trixi: CompressibleEulerEquations1D
using LinearAlgebra
using GZip
# used to generate the testing data for easy copy paste
function get_norms(arr)
    l1_norm = norm(arr, 1)
    l2_norm = norm(arr, 2)
    linf_norm = norm(arr, Inf)
    return l1_norm, l2_norm, linf_norm
end

function test_arr(arr, l1_norm, l2_norm, linf_norm; tol = 1e-14)
    @test isapprox(norm(arr, 1), l1_norm; atol = tol, rtol = tol)
    @test isapprox(norm(arr, 2), l2_norm; atol = tol, rtol = tol)
    @test isapprox(norm(arr, Inf), linf_norm; atol = tol, rtol = tol)
end

function set_lin_array!(arr::AbstractArray; normalizer::Real = 1.0)
    for I in CartesianIndices(arr)
        I_tuple = get_ijk(I) # I_tuple = (i, j, k)
        arr[I_tuple...] = (sum(I_tuple) + 10) * normalizer # arr[i, j, k] = (i + j + k)^2 + 10
    end
end

get_ijk(I::CartesianIndex{4}) = I[1], I[2], I[3], I[4]
get_ijk(I::CartesianIndex{3}) = I[1], I[2], I[3]
get_ijk(I::CartesianIndex{2}) = I[1], I[2]
get_ijk(I::CartesianIndex{1}) = I[1]

macro get_var_name(var)
    return QuoteNode(var)  # Returns the variable name as a Symbol
end

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
include("test_fluxes.jl")
include("test_update.jl")
