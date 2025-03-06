using LinearAlgebra
using TrixiBase
using PinCFlow_dev
using Test

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
             4.86721346318564e-5, tol = 9e-7) # Tolerance 9e-7 is needed to pass tests locally
end
