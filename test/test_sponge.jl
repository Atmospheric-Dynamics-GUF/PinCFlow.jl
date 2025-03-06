@testset "Test update.jl" begin
    semi = initialize_values(30, 1, 10, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    initialize_atmosphere!(semi)
    initialize_variables!(semi)
    (; cache, grid) = semi
    (; pStrat, rhoStrat, var, kr_sp_tfc, kr_sp_w_tfc, var0) = cache
    (; topography_surface) = grid

    normalizer = 0.001 # To get smaller values where arithmetic is more accurate
    set_lin_array_normalizer!(arr) = set_lin_array!(arr, normalizer = normalizer)
    set_lin_array_normalizer!.((pStrat, rhoStrat, topography_surface))

    set_lin_array_normalizer!(var.rho)
    set_lin_array_normalizer!(var.rhop)
    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)
    set_lin_array_normalizer!(var.exner)

    (; rho, rhop, u, v, w) = cache.var

    PinCFlow_dev.initialize_sponge!(semi)
    @show get_norms(semi.sponge.alpha)
    test_arr(semi.sponge.alpha, 100.65646473821992, 4.061806768248702, 0.2929262565775407,
             tol = 1e-14)

    PinCFlow_dev.applyUnifiedSponge_rho!(semi, dt)
    test_arr(rho, 148.7124917010242, 2.276604880034531, 0.06, tol = 1e-14)

    PinCFlow_dev.applyUnifiedSponge_rhop!(semi, dt)
    test_arr(rhop, 148.7124917010242, 2.276604880034531, 0.06, tol = 1e-14)

    PinCFlow_dev.applyUnifiedSponge_u!(semi, dt)
    test_arr(u, 153.42208281542247, 2.319471693928768, 0.06, tol = 1e-14)

    PinCFlow_dev.applyUnifiedSponge_v!(semi, dt)
    test_arr(v, 153.47600000000006, 2.320194859318576, 0.06, tol = 1e-14)

    PinCFlow_dev.applyUnifiedSponge_w!(semi, dt)
    test_arr(w, 153.47600000000006, 2.3201176568978767, 0.06, tol = 1e-14)
end
