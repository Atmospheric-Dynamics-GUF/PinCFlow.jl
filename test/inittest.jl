using PinCFlow_dev

semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)

initialize_atmosphere!(semi)

initialize_variables!(semi);

setBoundary!(semi)

semi.cache.var0.u .= semi.cache.var.u # TODO this is actually done in the time loop

compute_fluxes!(semi)
