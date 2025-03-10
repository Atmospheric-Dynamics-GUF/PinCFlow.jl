using TimerOutputs
using TrixiBase

function pincflow(model, dt)
    reset_timer!(timer())
    @trixi_timeit timer() "Pincflow" begin
    #! format: noindent
    initialize!(model)

    dt = dt / model.constants.tref

    ode = model.time

    setBoundary!(model)

    # TODO - These should be in parameters?
    errFlagBicg = false
    nIter = 0
    facprs = 1
    facray = 1

    Corrector(model, 1.0, errFlagBicg, nIter, "expl", facray, facprs)
    return
    initialize_sponge!(semi)

    time = 0.0
    for it in 1:100
        setBoundary!(semi)

        @trixi_timeit timer() "Misc" semi.cache.var.rhop .= semi.cache.var.rho
        @trixi_timeit timer() "Misc" copytuple_to_tuple!(semi.cache.var0,
            semi.cache.var)

        for RKStage in 1:3
            setBoundary!(semi)
            reconstruction!(semi)

            compute_fluxes!(semi)
            setBoundary_flux!(semi)

            @trixi_timeit timer() "Misc" if RKStage == 1
                copytuple_to_tuple_flux!(ode.flux0, semi.cache.flux)
            end

            @trixi_timeit timer() "Misc" ode.rhoOld .= semi.cache.var.rho

            massUpdate!(semi, ode, 0.5 * ode.dt, "rho", "tot", "expl", RKStage, 1)
            applyUnifiedSponge_rho!(semi, 0.5 * ode.stepFrac[RKStage] * dt)
            massUpdate!(semi, ode, 0.5 * ode.dt, "rhop", "lhs", "expl", RKStage, 1)
            applyUnifiedSponge_rhop!(semi, 0.5 * ode.stepFrac[RKStage] * dt)

            setBoundary!(semi)

            momentumPredictor!(semi, ode, 0.5 * dt, "lhs", "expl", RKStage, 1)

            applyUnifiedSponge_uvw!(semi, 0.5 * ode.stepFrac[RKStage] * dt)
        end
        RKStage = 4
        setBoundary!(semi)
        @trixi_timeit timer() "Misc" copytuple_to_tuple_flux!(semi.cache.flux,
            ode.flux0)

        @trixi_timeit timer() "Misc" ode.wOld .= semi.cache.var.w

        momentumPredictor!(semi, ode, 0.5 * dt, "rhs", "impl", RKStage, 1)

        setBoundary!(semi)

        massUpdate!(semi, ode, 0.5 * ode.dt, "rhop", "rhs", "impl", RKStage, 1)

        setBoundary!(semi)

        Corrector(semi, 0.5 * dt, errFlagBicg, nIter, "impl", 1, 1)

        setBoundary!(semi)

        @trixi_timeit timer() "Misc" begin
            copytuple_to_tuple!(ode.var1, semi.cache.var)
            copytuple_to_tuple_flux!(semi.cache.var, semi.cache.var0)
        end

        @trixi_timeit timer() "Misc" ode.rhopOld .= semi.cache.var.rhop

        massUpdate!(semi, ode, 0.5 * ode.dt, "rhop", "rhs", "expl", RKStage, 1)

        momentumPredictor!(semi, ode, 0.5 * dt, "rhs", "expl", RKStage, 1)

        setBoundary!(semi)

        copytuple_to_tuple!(semi.cache.var0, ode.var1)

        # TODO: ...
        ## set the boundary to var0

        for RKStage in 1:3
            setBoundary!(semi)

            reconstruction!(semi)

            #missing set boundary tilde var

            compute_fluxes!(semi)

            @trixi_timeit timer() "Boundary flux" setBoundary_flux!(semi)

            @trixi_timeit timer() "Misc" ode.rhoOld .= semi.cache.var.rho

            massUpdate!(semi, ode, ode.dt, "rho", "tot", "expl", RKStage, 1)

            applyUnifiedSponge_rho!(semi, ode.stepFrac[RKStage] * dt)

            massUpdate!(semi, ode, ode.dt, "rhop", "lhs", "expl", RKStage, 1)

            applyUnifiedSponge_rhop!(semi, ode.stepFrac[RKStage] * dt)

            setBoundary!(semi)

            momentumPredictor!(semi, ode, dt, "lhs", "expl", RKStage, 1)

            applyUnifiedSponge_uvw!(semi, ode.stepFrac[RKStage] * dt)
        end
        RKStage = 4
        setBoundary!(semi)
        @trixi_timeit timer() "Misc" copytuple_to_tuple_flux!(semi.cache.flux,
            ode.flux0)

        @trixi_timeit timer() "Misc" ode.wOld .= semi.cache.var.w

        momentumPredictor!(semi, ode, 0.5 * dt, "rhs", "impl", RKStage, 2)

        setBoundary!(semi)

        massUpdate!(semi, ode, 0.5 * ode.dt, "rhop", "rhs", "impl", RKStage, 2)

        setBoundary!(semi)

        Corrector(semi, 0.5 * dt, errFlagBicg, nIter, "impl", 2, 1)

        setBoundary!(semi)

        time = time + ode.dt
        @show time, it
    end
    end # timer
    display(timer())
end

function copytuple_to_tuple!(var, var1)
    var.u .= var1.u
    var.v .= var1.v
    var.w .= var1.w
    var.exner .= var1.exner
    var.rho .= var1.rho
    var.rhop .= var1.rhop
end

function copytuple_to_tuple_flux!(var, var1)
    var.u .= var1.u
    var.v .= var1.v
    var.w .= var1.w
    var.rho .= var1.rho
    var.rhop .= var1.rhop
end

function debugging!(semi)
    @show semi.cache.var.rho[5, 1, 5]
    @show semi.cache.var.rhop[5, 1, 5]
    @show semi.cache.var.u[5, 1, 5]
    @show semi.cache.var.v[5, 1, 5]
    @show semi.cache.var.w[5, 1, 5]
    @show semi.cache.var.exner[5, 1, 5]
    @show semi.cache.var.exner[5, 2, 1]
    @show semi.cache.var.u[5, 1, 1]
    @show semi.cache.var.u[5, 2, 1]
end
