function pincflow(semi, dt)

    initialize_atmosphere!(semi)

    dt = dt/semi.equations.tRef
    ode = time_discretization(semi, dt)

    initialize_variables!(semi);

    # add cleaning

    for RKstage = 1:3
      setBoundary!(semi)

      reconstruction!(semi)
 
      setBoundary!(semi)

      compute_fluxes!(semi)

      setBoundary_flux!(semi)

      if RKstage == 1
      #  ode.flux0 .= semi.cache.flux
      end

    # TODO: to be changed
      ode.rhoOld .= semi.cache.var.rho

      massUpdate_rho!(semi, ode, nothing, RKstage)

      ode.rhopOld .= semi.cache.var.rhop
      #call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, time, "rho")
      massUpdate_rhop!(semi, ode, "lhs", RKstage)

     # call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, time, "rhop")

      setBoundary!(semi)

      momentumPredictor!(semi, ode, "lhs", RKstage)

      #dt_Poisson = beta(RKstage) * dt

      #call Corrector(var, flux, dMom, dt_Poisson, errFlagBicg, nIterBicg, &
              #RKstage, "expl", 1., 1.)

     # call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, time, "uvw")
    end

    setBoundary!(semi)

    #semi.cache.flux .= semi.ode.flux0

    #semi.ode.wOld .= semi.cache.var.w

    #momentumPredictor!(semi, ode, "rhs", RKstage)



end