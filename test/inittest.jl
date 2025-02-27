using PinCFlow_dev

semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)

initialize_atmosphere!(semi)

initialize_variables!(semi);

setBoundary!(semi)

semi.cache.var0.u .= semi.cache.var.u # TODO this is actually done in the time loop

compute_fluxes!(semi)

dt = 1.0
ode = time_discretization(semi, dt)

massUpdate!(semi, ode, "lhs", 2)
massUpdate!(semi, ode, "rhs", 1)

momentumPredictor_u!(semi, ode, "lhs", 1)
momentumPredictor_u!(semi, ode, "rhs", 1)

momentumPredictor_v!(semi, ode, "lhs", 1)
momentumPredictor_v!(semi, ode, "rhs", 1)

momentumPredictor_w!(semi, ode, "lhs", 1)
momentumPredictor_w!(semi, ode, "rhs", 1)