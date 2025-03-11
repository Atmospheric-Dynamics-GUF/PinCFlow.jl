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
    display(timer())
    initialize_sponge!(model)


    time = 0.0
    for it in 1:100
        setBoundary!(model)
        @trixi_timeit timer() "Misc" model.variables.prognostic_fields.rhop .= model.variables.prognostic_fields.rho
        @trixi_timeit timer() "Misc" copytuple_to_tuple!(model.variables.prognostic_fields_0,
            model.variables.prognostic_fields)
        flux0 = copy(model.fluxes)
        for RKStage in 1:3
            setBoundary!(model)
            reconstruction!(model)

            compute_fluxes!(model)
            setBoundary_flux!(model)

            @trixi_timeit timer() "Misc" if RKStage == 1
                flux0 = copy(model.fluxes)
                # copytuple_to_tuple_flux!(ode.flux0, semi.cache.flux)
            end

            @trixi_timeit timer() "Misc" model.variables.history.rho .= model.variables.prognostic_fields.rho

            massUpdate!(model, ode, 0.5 * dt, "rho", "tot", "expl", RKStage, 1)
            applyUnifiedSponge_rho!(model, 0.5 * ode.stepfrac[RKStage] * dt)
            massUpdate!(model, ode, 0.5 * dt, "rhop", "lhs", "expl", RKStage, 1)
            applyUnifiedSponge_rhop!(model, 0.5 * ode.stepfrac[RKStage] * dt)
            setBoundary!(model)

            momentumPredictor!(model, ode, 0.5 * dt, "lhs", "expl", RKStage, 1)
            applyUnifiedSponge_uvw!(model, 0.5 * ode.stepfrac[RKStage] * dt)
        end
        RKStage = 4
        setBoundary!(model)
        @trixi_timeit timer() "Misc" model.fluxes = copy(flux0)
        # @trixi_timeit timer() "Misc" copytuple_to_tuple_flux!(semi.cache.flux,
        #     ode.flux0)

        @trixi_timeit timer() "Misc" model.variables.history.w .= model.variables.prognostic_fields.w

        momentumPredictor!(model, ode, 0.5 * dt, "rhs", "impl", RKStage, 1)

        setBoundary!(model)

        massUpdate!(model, ode, 0.5 * dt, "rhop", "rhs", "impl", RKStage, 1)

        setBoundary!(model)

        Corrector(model, 0.5 * dt, errFlagBicg, nIter, "impl", 1, 1)

        setBoundary!(model)

        @trixi_timeit timer() "Misc" begin
            copytuple_to_tuple!(model.variables.prognostic_fields_1, model.variables.prognostic_fields)
            copytuple_to_tuple!(model.variables.prognostic_fields, model.variables.prognostic_fields_0)
        end

        @trixi_timeit timer() "Misc" model.variables.history.rhop .= model.variables.prognostic_fields.rhop

        massUpdate!(model, ode, 0.5 * dt, "rhop", "rhs", "expl", RKStage, 1)

        momentumPredictor!(model, ode, 0.5 * dt, "rhs", "expl", RKStage, 1)

        setBoundary!(model)

        copytuple_to_tuple!(model.variables.prognostic_fields_0, model.variables.prognostic_fields_1)

        # TODO: ...
        ## set the boundary to var0

        for RKStage in 1:3
            setBoundary!(model)

            reconstruction!(model)

            #missing set boundary tilde var

            compute_fluxes!(model)

            @trixi_timeit timer() "Boundary flux" setBoundary_flux!(model)

            @trixi_timeit timer() "Misc" model.variables.history.rho .= model.variables.prognostic_fields.rho

            massUpdate!(model, ode, dt, "rho", "tot", "expl", RKStage, 1)

            applyUnifiedSponge_rho!(model, ode.stepfrac[RKStage] * dt)

            massUpdate!(model, ode, dt, "rhop", "lhs", "expl", RKStage, 1)

            applyUnifiedSponge_rhop!(model, ode.stepfrac[RKStage] * dt)

            # TODO ode.dt == dt?
            setBoundary!(model)

            momentumPredictor!(model, ode, dt, "lhs", "expl", RKStage, 1)

            applyUnifiedSponge_uvw!(model, ode.stepfrac[RKStage] * dt)
        end
        RKStage = 4
        setBoundary!(model)
        @trixi_timeit timer() "Misc" model.fluxes = copy(flux0)

        @trixi_timeit timer() "Misc" model.variables.history.w .= model.variables.prognostic_fields.w

        momentumPredictor!(model, ode, 0.5 * dt, "rhs", "impl", RKStage, 2)

        setBoundary!(model)

        massUpdate!(model, ode, 0.5 * dt, "rhop", "rhs", "impl", RKStage, 2)

        setBoundary!(model)

        Corrector(model, 0.5 * dt, errFlagBicg, nIter, "impl", 2, 1)

        setBoundary!(model)

        time = time + dt
        @show time, it
    end
    end # timer
    display(timer())
end

function copytuple_to_tuple!(var, var1)
    var.u .= var1.u
    var.v .= var1.v
    var.w .= var1.w
    var.pip .= var1.pip
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
