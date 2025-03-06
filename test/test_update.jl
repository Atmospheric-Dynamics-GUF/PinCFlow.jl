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
    ode = time_discretization(semi, 0.1)

    set_lin_array_normalizer!(ode.dRho)
    set_lin_array_normalizer!(ode.dRhop)
    set_lin_array_normalizer!(ode.dMom)

    (; rho, rhop, u, v, w) = cache.var

    PinCFlow_dev.massUpdate!(semi, ode, 0.5 * 1.0, "rho", "tot", "expl", 2, 1)
    test_arr(rho, 147.2038515629778, 2.2610408292790254, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(ode.dRho)

    set_lin_array_normalizer!(rho)
    PinCFlow_dev.massUpdate!(semi, ode, 0.5 * 1.0, "rhop", "rhs", "impl", 2, 1)
    test_arr(rhop, 153.18651104551827, 2.319981249089708, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(ode.dRhop)

    set_lin_array_normalizer!(rhop)
    PinCFlow_dev.massUpdate!(semi, ode, 0.5 * 1.0, "rhop", "lhs", "expl", 2, 1)
    test_arr(rhop, 147.2038515629778, 2.2610408292790254, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(ode.dRhop)

    set_lin_array_normalizer!(rhop)
    PinCFlow_dev.massUpdate!(semi, ode, 0.5 * 1.0, "rhop", "rhs", "expl", 2, 1)
    test_arr(rhop, 153.52198373074026, 2.3249185833765122, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(rhop)
    set_lin_array_normalizer!(ode.dRhop)

    set_lin_array_normalizer!(u)
    set_lin_array_normalizer!(v)
    set_lin_array_normalizer!(w)
    set_lin_array_normalizer!(ode.usave)
    PinCFlow_dev.momentumPredictor!(semi, ode, 0.5 * 1.0, "lhs", "expl", 2, 1)
    +
    test_arr(u, 238.74929329766223, 5.881635566827642, 0.4273353001174364, tol = 1e-14)
    test_arr(v, 323.59998405004694, 8.07055523588239, 0.44816863345076985, tol = 1e-14)
    test_arr(w, 232.5328808987527, 5.824225035198435, 0.451871006844049, tol = 1e-14)
    test_arr(ode.usave, 153.47600000000006, 2.3241884605169103, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(u)
    set_lin_array_normalizer!(v)
    set_lin_array_normalizer!(w)
    set_lin_array_normalizer!(ode.usave)
    set_lin_array_normalizer!(ode.dMom)

    PinCFlow_dev.momentumPredictor!(semi, ode, 0.5 * 1.0, "rhs", "expl", 2, 1)

    test_arr(u, 152.29547184017403, 2.309105007878862, 0.06, tol = 1e-14)
    test_arr(v, 153.36175533913394, 2.3226440800408104, 0.06, tol = 1e-14)
    test_arr(w, 205.80882774255787, 4.340259513308493, 0.24215723191622596, tol = 1e-14)
    test_arr(ode.usave, 152.29547184017403, 2.309105007878862, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(u)
    set_lin_array_normalizer!(v)
    set_lin_array_normalizer!(w)
    set_lin_array_normalizer!(ode.usave)
    set_lin_array_normalizer!(ode.dMom)
    PinCFlow_dev.momentumPredictor!(semi, ode, 0.5 * 1.0, "rhs", "impl", 2, 1)
    @show get_norms(u)
    @show get_norms(v)
    @show get_norms(w)
    @show get_norms(ode.usave)
    test_arr(u, 152.29547184017403, 2.309105007878862, 0.06, tol = 1e-14)
    test_arr(v, 153.36175533913394, 2.3226440800408104, 0.06, tol = 1e-14)
    test_arr(w, 203.72238988780234, 4.232319767715575, 0.2338514844262595, tol = 1e-14)
    test_arr(ode.usave, 152.29547184017403, 2.309105007878862, 0.06, tol = 1e-14)
end
