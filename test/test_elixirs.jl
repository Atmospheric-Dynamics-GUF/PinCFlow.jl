using LinearAlgebra
using TrixiBase
using PinCFlow_dev
using Test

# used to generate the testing data for easy copy paste
function get_norms(arr)
    l1_norm = norm(arr, 1)
    l2_norm = norm(arr, 2)
    linf_norm = norm(arr, Inf)
    return l1_norm, l2_norm, linf_norm
end

function test_arr(arr, l1_norm, l2_norm, linf_norm; atol = 1e-14, rtol = 1e-14)
    @test isapprox(norm(arr, 1), l1_norm; atol = atol, rtol = rtol)
    @test isapprox(norm(arr, 2), l2_norm; atol = atol, rtol = rtol)
    @test isapprox(norm(arr, Inf), linf_norm; atol = atol, rtol = rtol)
end

@testset "elixir_mountainwave" begin
    trixi_include(joinpath(pincflow_examples_dir(), "elixir_mountainwave.jl"),
                  nx = 9, nz = 3)

    test_arr(semi.cache.var.rho, 0.2932986276870159, 0.016954799596920693,
             0.0021413112274175743)
    test_arr(semi.cache.var.u, 34.36787057152193, 1.0829898397682243, 0.03657494978683221)
    test_arr(semi.cache.var.v, 0.0, 0.0, 0.0)
    test_arr(semi.cache.var.w, 0.24906743875884985, 0.014125193256054751,
             0.0014552046070389162)
    test_arr(semi.cache.var.rhop, 0.2626779252704753, 0.015356203279767772,
             0.0019438768803848995)
    test_arr(semi.cache.var.exner, 0.004415532955941557, 0.0003705079602120867,
             4.86721346318564e-5) # Only passes on github CI. Lower tolerance to 9e-7 to pass locally
end
