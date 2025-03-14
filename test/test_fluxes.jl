@testset "Test fluxes.jl" begin
    semi = initialize_values(30, 1, 10, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    initialize_atmosphere!(semi)
    initialize_variables!(semi)
    (; cache, grid) = semi
    rhs = cache.rhs_bicg
    (; pStrat, rhoStrat, var, kr_sp_tfc, kr_sp_w_tfc, var0) = cache
    (; topography_surface) = grid

    normalizer = 0.01 # To get smaller values where arithmetic is more accurate
    set_lin_array_normalizer!(arr) = set_lin_array!(arr, normalizer = normalizer)
    set_lin_array_normalizer!.((pStrat, rhoStrat, topography_surface))

    set_lin_array_normalizer!(var.rho)
    set_lin_array_normalizer!(var.rhop)
    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)

    PinCFlow.reconstruction!(semi)
    (; rhoTilde, rhopTilde, uTilde, vTilde, wTilde) = cache

    test_arr(rhoTilde, 15120.0, 122.96340919151518, 1.0, tol = 1e-14)
    test_arr(rhopTilde, 15120.0, 122.96340919151518, 1.0, tol = 1e-14)
    test_arr(uTilde, 9374.400000000001, 80.74848357709314, 1.1300000000000001, tol = 1e-14)
    test_arr(vTilde, 9374.400000000001, 80.74805260809762, 1.1300000000000001, tol = 1e-14)
    test_arr(wTilde, 8703.234242111952, 73.30783802039744, 0.8492267931815118, tol = 1e-14)

    set_lin_array_normalizer!(var0.rho)
    set_lin_array_normalizer!(var0.rhop)
    set_lin_array_normalizer!(var0.u)
    set_lin_array_normalizer!(var0.v)
    set_lin_array_normalizer!(var0.w)
    PinCFlow.compute_fluxes!(semi)
    (; flux) = cache

    test_arr(flux.rho, 234.48977709288997, 7.465169904723465, 0.4313585705042774,
             tol = 1e-14)
    test_arr(flux.rhop, 117.24488854644498, 3.7325849523617327, 0.2156792852521387,
             tol = 1e-14)
    test_arr(flux.u, 88.4864783397629, 3.0875380464248465, 0.22591343710336997, tol = 1e-14)
    test_arr(flux.v, 149.20376149456763, 3.964286453926213, 0.22591343710336997,
             tol = 1e-14)
    test_arr(flux.w, 81.05895692784512, 2.705278309488443, 0.18549395063610127, tol = 1e-14)
end
