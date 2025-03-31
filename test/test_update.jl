using Setfield

@testset "Test update.jl" begin
    p = default_parameters()
    domain = DomainParameters(sizex = 30, sizez = 10, nbx = 3, nby = 3, nbz = 3)
    p = @set p.domain = domain

    model = Model(p)
    initialize!(model)

    var = model.variables.prognostic_fields
    flux = model.fluxes
    topography_surface = model.grid.topography_surface
    rhs = model.operator.cache.rhs_bigc
    (; pstrattfc, rhostrattfc) = model.atmosphere
    (; kr_sp_tfc, kr_sp_w_tfc) = model.atmosphere.sponge
    var0 = model.variables.prognostic_fields_0
    ode = model.time
    dt = 0.1

    normalizer = 0.001 # To get smaller values where arithmetic is more accurate
    set_lin_array_normalizer!(arr) = set_lin_array!(arr, normalizer = normalizer)
    set_lin_array_normalizer!.((pstrattfc, rhostrattfc, topography_surface))

    set_lin_array_normalizer!(var.rho)
    set_lin_array_normalizer!(var.rhop)
    set_lin_array_normalizer!(var.u)
    set_lin_array_normalizer!(var.v)
    set_lin_array_normalizer!(var.w)
    set_lin_array_normalizer!(var.pip)

    set_lin_array_normalizer!(flux.rho)
    set_lin_array_normalizer!(flux.rhop)
    set_lin_array_normalizer!(flux.u)
    set_lin_array_normalizer!(flux.v)
    set_lin_array_normalizer!(flux.w)

    #TODO: in the old code, ode = time_discretization is called here. this sets uOld, vOld, wOld
    model.variables.history.u .= var.u
    model.variables.history.v .= var.v
    model.variables.history.w .= var.w
    # this is obviously wrong. leave for now to pass the tests
    model.variables.history.rho .= var.u
    model.variables.history.rhop .= var.u

    set_lin_array_normalizer!(model.variables.tendencies.drho)
    set_lin_array_normalizer!(model.variables.tendencies.drhop)
    set_lin_array_normalizer!(model.variables.tendencies.dmom)

    (; rho, rhop, u, v, w) = var

    PinCFlow.massUpdate!(model, ode, 0.5 * 1.0, "rho", "tot", "expl", 2, 1)
    test_arr(rho, 147.2038515629778, 2.2610408292790254, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(model.variables.tendencies.drho)

    set_lin_array_normalizer!(rho)
    PinCFlow.massUpdate!(model, ode, 0.5 * 1.0, "rhop", "rhs", "impl", 2, 1)
    test_arr(rhop, 153.18651104551827, 2.319981249089708, 0.06, tol = 1e-14)

    set_lin_array_normalizer!(model.variables.tendencies.drhop)

    set_lin_array_normalizer!(rhop)
    PinCFlow.massUpdate!(model, ode, 0.5 * 1.0, "rhop", "lhs", "expl", 2, 1)
    test_arr(rhop, 147.2038515629778, 2.2610408292790254, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(model.variables.tendencies.drhop)

    set_lin_array_normalizer!(rhop)
    PinCFlow.massUpdate!(model, ode, 0.5 * 1.0, "rhop", "rhs", "expl", 2, 1)
    test_arr(rhop, 153.52198373074026, 2.3249185833765122, 0.06, tol = 1e-14)
    set_lin_array_normalizer!(rhop)
    set_lin_array_normalizer!(model.variables.tendencies.drhop)

    set_lin_array_normalizer!(u)
    set_lin_array_normalizer!(v)
    set_lin_array_normalizer!(w)
    set_lin_array_normalizer!(model.variables.usave)

    PinCFlow.momentumPredictor!(model, ode, 0.5 * 1.0, "lhs", "expl", 2, 1)

    test_arr(u, 238.74929329766223, 5.881635566827642, 0.4273353001174364, tol = 1e-14)
    test_arr(v, 323.59998405004694, 8.07055523588239, 0.44816863345076985, tol = 1e-14)
    test_arr(w, 232.5328808987527, 5.824225035198435, 0.451871006844049, tol = 1e-14)
    test_arr(model.variables.usave, 153.47600000000006, 2.3241884605169103, 0.06,
             tol = 1e-14)

    set_lin_array_normalizer!(u)
    set_lin_array_normalizer!(v)
    set_lin_array_normalizer!(w)
    set_lin_array_normalizer!(model.variables.usave)
    set_lin_array_normalizer!(model.variables.tendencies.dmom)

    PinCFlow.momentumPredictor!(model, ode, 0.5 * 1.0, "rhs", "expl", 2, 1)

    test_arr(u, 152.29547184017403, 2.309105007878862, 0.06, tol = 1e-14)
    test_arr(v, 153.36175533913394, 2.3226440800408104, 0.06, tol = 1e-14)
    test_arr(w, 205.80882774255787, 4.340259513308493, 0.24215723191622596, tol = 1e-14)
    test_arr(model.variables.usave, 152.29547184017403, 2.309105007878862, 0.06,
             tol = 1e-14)

    set_lin_array_normalizer!(u)
    set_lin_array_normalizer!(v)
    set_lin_array_normalizer!(w)
    set_lin_array_normalizer!(model.variables.usave)
    set_lin_array_normalizer!(model.variables.tendencies.dmom)
    PinCFlow.momentumPredictor!(model, ode, 0.5 * 1.0, "rhs", "impl", 2, 1)

    test_arr(u, 152.29547184017403, 2.309105007878862, 0.06, tol = 1e-14)
    test_arr(v, 153.36175533913394, 2.3226440800408104, 0.06, tol = 1e-14)
    test_arr(w, 203.72238988780234, 4.232319767715575, 0.2338514844262595, tol = 1e-14)
    test_arr(model.variables.usave, 152.29547184017403, 2.309105007878862, 0.06,
             tol = 1e-14)
end
