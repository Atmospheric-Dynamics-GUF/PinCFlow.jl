using LinearAlgebra
using TrixiBase
using PinCFlow
using Test

@testset "elixir_mountainwave" begin
    trixi_include(joinpath(pincflow_examples_dir(), "elixir_mountainwave.jl"),
                  nx = 9, nz = 3)

    test_arr(semi.cache.var.rho, 0.30291020305314476, 0.018067350880649553,
             0.002329724957484748)
    test_arr(semi.cache.var.u, 34.798656811192025, 1.0968034985215316, 0.03709804969480197)
    test_arr(semi.cache.var.v, 0.0, 0.0, 0.0)
    test_arr(semi.cache.var.w, 0.21756140137542518, 0.012481925847136642,
             0.0013221551824446711)
    test_arr(semi.cache.var.rhop, 0.27473000166815187, 0.016476285614253763,
             0.002122855224708782)
    test_arr(semi.cache.var.exner, 0.018105619057978757, 0.001428064135088831,
             0.00013806311680403024, tol = 9e-7) # Tolerance 9e-7 is needed to pass tests locally
end
