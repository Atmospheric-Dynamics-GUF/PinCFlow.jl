@testset "Test update.jl" begin
    semi = initialize_values(30, 1, 10, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    initialize_atmosphere!(semi)
    initialize_variables!(semi)
    (; cache, grid) = semi
    rhs = cache.rhs_bicg
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

    set_lin_array_normalizer!(cache.flux.rho)
    set_lin_array_normalizer!(cache.flux.rhop)
    set_lin_array_normalizer!(cache.flux.u)
    set_lin_array_normalizer!(cache.flux.v)
    set_lin_array_normalizer!(cache.flux.w)

    (; rho, rhop, u, v, w, exner) = cache.var

    PinCFlow_dev.setBoundary_x!(semi, PinCFlow_dev.PeriodicBC())

    test_arr(u, 157.55600000000004, 2.3564405360628085, 0.057, tol = 1e-14)
    test_arr(v, 153.47600000000006, 2.3083353309257246, 0.057, tol = 1e-14)
    test_arr(w, 153.47600000000006, 2.3083353309257246, 0.057, tol = 1e-14)
    test_arr(rho, 153.47600000000006, 2.3083353309257246, 0.057, tol = 1e-14)
    test_arr(rhop, 153.47600000000006, 2.3083353309257246, 0.057, tol = 1e-14)
    test_arr(exner, 153.47600000000003, 2.3224323456238647, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(var.rho)
    set_lin_array_normalizer!(var.rhop)
    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)
    set_lin_array_normalizer!(var.exner)
    PinCFlow_dev.setBoundary_y!(semi, PinCFlow_dev.PeriodicBC())

    test_arr(u, 153.47600000000006, 2.3203965178391304, 0.057, tol = 1e-14)
    test_arr(v, 155.99200000000005, 2.3516309234231376, 0.057, tol = 1e-14)
    test_arr(w, 153.47600000000006, 2.3203965178391304, 0.057, tol = 1e-14)
    test_arr(rho, 153.47600000000006, 2.3203965178391304, 0.057, tol = 1e-14)
    test_arr(rhop, 153.47600000000006, 2.3203965178391304, 0.057, tol = 1e-14)
    test_arr(exner, 153.47600000000006, 2.3239178126603286, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(var.rho)
    set_lin_array_normalizer!(var.rhop)
    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)
    set_lin_array_normalizer!(var.exner)
    PinCFlow_dev.setBoundary_z!(semi, PinCFlow_dev.PeriodicBC())

    test_arr(u, 153.47600000000006, 2.3126979915241814, 0.057, tol = 1e-14)
    test_arr(v, 153.47600000000006, 2.3126979915241814, 0.057, tol = 1e-14)
    test_arr(w, 156.43600000000004, 2.3476311464963966, 0.057, tol = 1e-14)
    test_arr(rho, 153.47600000000006, 2.3126979915241814, 0.057, tol = 1e-14)
    test_arr(rhop, 153.47600000000006, 2.3126979915241814, 0.057, tol = 1e-14)
    test_arr(exner, 153.47600000000006, 2.3229145485789973, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(var.rho)
    set_lin_array_normalizer!(var.rhop)
    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)
    set_lin_array_normalizer!(var.exner)
    PinCFlow_dev.setBoundary_z!(semi, PinCFlow_dev.SolidWallBC())

    test_arr(u, 153.47600000000006, 2.31269799152418, 0.057, tol = 1e-14)
    test_arr(v, 153.47600000000006, 2.31269799152418, 0.057, tol = 1e-14)
    test_arr(w, 135.42000000000004, 2.166727486325857, 0.056, tol = 1e-14)
    test_arr(rho, 153.47600000000006, 2.31269799152418, 0.057, tol = 1e-14)
    test_arr(rhop, 153.47600000000006, 2.31269799152418, 0.057, tol = 1e-14)
    test_arr(exner, 153.47600000000006, 2.3229145485789977, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(cache.flux.w)
    PinCFlow_dev.setBoundary_flux!(semi)
end
