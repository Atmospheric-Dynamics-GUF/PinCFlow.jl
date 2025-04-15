function integrate(namelists::Namelists)

    #-------------------------------------------------
    #                    Set up
    #-------------------------------------------------

    # Initialize time variables.
    itime = 0
    rkstage = 0
    time = 0.0
    dt = 0.0

    # Initialize I/O variables.
    iout = 0
    output = false
    nextoutputtime = 0.0

    # Initialize Poisson variables.
    errflagbicg = false
    niterbicg = 0
    ntotalbicg = 0
    naveragebicg = 0.0

    # Initialize the model state.
    state = State(namelists)

    # Save CPU start time.
    cpu_start_time = now()

    # Get all necessary fields.
    (; nprocx, nprocy) = state.namelists.domain
    (; model) = state.namelists.setting
    (; initialcleaning) = state.namelists.poisson
    (; dtmin_dim) = state.namelists.discretization
    (; restart, maxtime, outputtimediff, output_steps, maxiter, noutput) =
        state.namelists.output
    (; nstages, stepfrac) = state.time
    (; tref) = state.constants
    (; master) = state.domain

    # Print information.
    if master
        println("")
        println(repeat("-", 80))
        println(repeat(" ", 36), "PinCFlow")
        println(
            repeat(" ", 12),
            "developed by Rieper et al (2013) and Schmid et al (2021)",
        )
        println(repeat(" ", 28), "modified by many others")
        println(repeat("-", 80))
        println("")
        println("Date: ", Dates.Date(cpu_start_time))
        println("Time: ", Dates.Time(cpu_start_time))
        println("")
        println("Virtual topology: (nprocx, nprocy) = ", (nprocx, nprocy))
        println("")
    end

    #---------------------------------------------
    #        Initial divergence cleaning
    #---------------------------------------------

    if initialcleaning
        set_boundaries!(state, BoundaryPredictands())

        (errflagbicg, niterbicg) =
            apply_corrector!(state, 1.0, EXPL(), 1.0, 1.0)

        if errflagbicg
            exit()
        end
    end

    #---------------------------------------------
    #              Initialize MS-GWaM
    #---------------------------------------------

    initialize_rays!(state)

    #-------------------------------------------------
    #              Read initial data
    #-------------------------------------------------

    if restart
        if master
            println("Reading restart file...")
        end

        time = read_input!(state)

        if maxtime < time * tref
            println("Restart error: maxtime < time!")
        end

        set_boundaries!(state, BoundaryPredictands())
    end

    #------------------------------------------
    #              Initial output
    #------------------------------------------

    # Create the output file.
    create_output(state)

    # Write the initial state.
    iout = write_output(state, time, iout, cpu_start_time)

    output = false
    nextoutputtime = time * tref + outputtimediff

    #-----------------------------------------------------
    #                        Time loop
    #-----------------------------------------------------

    if master
        println("Starting the time loop...")
    end

    if output_steps
        maxiterations = maxiter
    else
        maxiterations = 2^30
    end

    for itime in 1:maxiterations
        if master
            println("")
            println(repeat("-", 80))
            println("Time step = ", itime)
            println("Time = ", time * tref, " seconds")
            println(repeat("-", 80))
        end

        #----------------------------------
        #         Calc time step
        #----------------------------------

        dt = compute_time_step(state)

        # correct dt to hit desired output time for outputType 'time'
        if !output_steps
            if (time + dt) * tref + dtmin_dim > nextoutputtime
                dt = nextoutputtime / tref - time
                output = true
                if master
                    println(
                        "Time step for output: dt = ",
                        dt * tref,
                        " seconds",
                    )
                end
            end
        end

        # Advance time
        time += dt

        #-----------------------------------------------------------------
        #                         Sponge layer
        #-----------------------------------------------------------------

        compute_sponge!(state, dt)

        #-----------------------------------------------------------------
        #                       Apply saturation scheme
        #-----------------------------------------------------------------

        apply_saturation_scheme!(state, dt)

        #---------------------------------------------------------------
        #                   Semi-implicit time scheme
        #---------------------------------------------------------------

        # Run MS-GWaM.
        for rkstage in 1:nstages
            propagate_rays!(state, dt, rkstage)
        end
        split_rays!(state)
        shift_rays!(state)
        merge_rays!(state)
        set_boundary_rays!(state)
        compute_mean_flow_effect!(state)

        # Synchronization of density fluctuations
        synchronize_density_fluctuations!(state)

        set_boundaries!(state, BoundaryPredictands())

        # put initial state into var0 in order to save the advecting
        # velocities
        v0 = deepcopy(state.variables.predictands)
        f0 = deepcopy(state.variables.fluxes)

        # (1) explicit integration of convective and
        #     viscous-diffusive/turbulent fluxes over half a time step,
        #     with the advection velocity kept constant
        #     \psi^# = \psi^n + A^{dt/2} (\psi^n, v^n)

        if master
            println("Beginning a semi-implicit time step...")
            println("(1) Explicit integration lhs over dt/2...")
        end

        for rkstage in 1:nstages
            set_boundaries!(state, BoundaryPredictands())

            # Reconstruction
            reconstruct!(state)
            set_boundaries!(state, BoundaryReconstructions())

            # Fluxes
            compute_fluxes!(state, v0)
            set_boundaries!(state, BoundaryFluxes())

            # store initial flux
            if rkstage == 1
                f0 = deepcopy(state.variables.fluxes)
            end

            # RK step for density and density fluctuations

            state.variables.backups.rhoold .= state.variables.predictands.rho

            update!(state, 0.5 * dt, rkstage, Rho(), LHS(), EXPL())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                Rho(),
                model,
            )

            update!(state, 0.5 * dt, rkstage, RhoP(), LHS(), EXPL())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                RhoP(),
                model,
            )

            # RK step for momentum

            set_boundaries!(state, BoundaryPredictands())

            update!(state, 0.5 * dt, rkstage, U(), LHS(), EXPL())
            update!(state, 0.5 * dt, rkstage, V(), LHS(), EXPL())
            update!(state, 0.5 * dt, rkstage, W(), LHS(), EXPL())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                U(),
                model,
            )
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                V(),
                model,
            )
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                W(),
                model,
            )
        end

        # (2) implicit integration of the linear right-hand sides of the
        #     equations for density fluctuations and momentum over half a
        #     time step, under consideration of the divergence constraint
        #     \psi^{n+1/2} = \psi^# + dt/2 Q(\psi^{n+1/2})

        if master
            println("(2) Implicit integration rhs over dt/2...")
        end

        set_boundaries!(state, BoundaryPredictands())

        # use initial flux for update of reference atmosphere and w0
        state.variables.fluxes.phirho .= f0.phirho
        state.variables.fluxes.phirhop .= f0.phirhop
        state.variables.fluxes.phiu .= f0.phiu
        state.variables.fluxes.phiv .= f0.phiv
        state.variables.fluxes.phiw .= f0.phiw

        # uStar and vStar are needed for update of density fluctuations,
        # therefore w is stored instead of rhop

        state.variables.backups.wold .= state.variables.predictands.w

        # update winds (uStar, vStar, wStar)

        update!(state, 0.5 * dt, U(), RHS(), IMPL(), 1.0)
        update!(state, 0.5 * dt, V(), RHS(), IMPL(), 1.0)
        update!(state, 0.5 * dt, W(), RHS(), IMPL(), 1.0)

        set_boundaries!(state, BoundaryPredictands())

        # update density fluctuations (rhopStar)

        update!(state, 0.5 * dt, RhoP(), RHS(), IMPL(), 1.0)

        set_boundaries!(state, BoundaryPredictands())

        # Correct momentum and density fluctuations
        (errflagbicg, niterbicg) =
            apply_corrector!(state, 0.5 * dt, IMPL(), 1.0, 1.0)

        if errflagbicg
            iout = write_output(state, time, iout, cpu_start_time)
            if master
                println("Output last state into record", iout)
            end
            exit()
        end

        ntotalbicg += niterbicg

        set_boundaries!(state, BoundaryPredictands())

        # put new state into var1 in order to save the advection velocities
        v1 = deepcopy(state.variables.predictands)

        # (3) explicit integration of the linear right-hand sides of the
        #     equations for density fluctuations and momentum over half a
        #     time step, under consideration of the divergence constraint
        #     \psi^\ast = \psi^n + dt/2 Q(\psi^n)

        if master
            println("(3) Explicit integration rhs over dt/2...")
        end

        # (3) uses updated pressure field and (5) adjusts pressure over half a
        # time step#
        state.variables.predictands.rho .= v0.rho
        state.variables.predictands.rhop .= v0.rhop
        state.variables.predictands.u .= v0.u
        state.variables.predictands.v .= v0.v
        state.variables.predictands.w .= v0.w

        set_boundaries!(state, BoundaryPredictands())

        state.variables.backups.rhopold .= state.variables.predictands.rhop

        # update density fluctuations (rhopStar)
        update!(state, 0.5 * dt, RhoP(), RHS(), EXPL())

        # update winds (uStar, vStar, wStar)

        update!(state, 0.5 * dt, U(), RHS(), EXPL())
        update!(state, 0.5 * dt, V(), RHS(), EXPL())
        update!(state, 0.5 * dt, W(), RHS(), EXPL())

        set_boundaries!(state, BoundaryPredictands())

        # (4) explicit integration of convective and
        #     viscous-diffusive/turbulent fluxes over a full time step,
        #     with the advection velocity kept constant
        #     \psi^{\ast\ast} = \psi^\ast + A^dt (\psi^\ast, v^{n+1/2})

        if master
            println("(4) Explicit integration lhs over dt...")
        end

        v0 = deepcopy(v1)

        for rkstage in 1:nstages
            set_boundaries!(state, BoundaryPredictands())

            # Reconstruction
            reconstruct!(state)
            set_boundaries!(state, BoundaryReconstructions())

            # Fluxes
            compute_fluxes!(state, v0)
            set_boundaries!(state, BoundaryFluxes())

            # RK step for density and density fluctuations

            state.variables.backups.rhoold .= state.variables.predictands.rho

            update!(state, dt, rkstage, Rho(), LHS(), EXPL())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * dt,
                time,
                Rho(),
                model,
            )

            update!(state, dt, rkstage, RhoP(), LHS(), EXPL())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * dt,
                time,
                RhoP(),
                model,
            )

            # RK step for momentum
            set_boundaries!(state, BoundaryPredictands())

            update!(state, dt, rkstage, U(), LHS(), EXPL())
            update!(state, dt, rkstage, V(), LHS(), EXPL())
            update!(state, dt, rkstage, W(), LHS(), EXPL())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * dt,
                time,
                U(),
                model,
            )
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * dt,
                time,
                V(),
                model,
            )
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * dt,
                time,
                W(),
                model,
            )
        end

        # (5) implicit integration of the linear right-hand sides of the
        #     equations for density fluctuations and momentum over half a
        #     time step, under consideration of the divergence constraint
        #     \psi^{n+1} = \psi^{\ast\ast} + dt/2 Q(\psi^{n+1})

        if master
            println("(5) Implicit integration rhs over dt/2...")
        end

        set_boundaries!(state, BoundaryPredictands())

        # use initial flux for update of reference atmosphere and w0
        state.variables.fluxes.phirho .= f0.phirho
        state.variables.fluxes.phirhop .= f0.phirhop
        state.variables.fluxes.phiu .= f0.phiu
        state.variables.fluxes.phiv .= f0.phiv
        state.variables.fluxes.phiw .= f0.phiw

        # uStar and vStar are needed for update of density fluctuations,
        # therefore w is stored instead of rhop
        state.variables.backups.wold .= state.variables.predictands.w

        # update winds (uStar, vStar, wStar)

        update!(state, 0.5 * dt, U(), RHS(), IMPL(), 2.0)
        update!(state, 0.5 * dt, V(), RHS(), IMPL(), 2.0)
        update!(state, 0.5 * dt, W(), RHS(), IMPL(), 2.0)

        set_boundaries!(state, BoundaryPredictands())

        # update density fluctuations (rhopStar)

        update!(state, 0.5 * dt, RhoP(), RHS(), IMPL(), 2.0)

        set_boundaries!(state, BoundaryPredictands())

        # (3) uses updated pressure field and (5) adjusts pressure over half a
        # time step
        (errflagbicg, niterbicg) =
            apply_corrector!(state, 0.5 * dt, IMPL(), 2.0, 1.0)

        if errflagbicg
            iout = write_output(state, time, iout, cpu_start_time)
            if master
                println("Output last state into record ", iout)
            end
            exit()
        end

        set_boundaries!(state, BoundaryPredictands())

        ntotalbicg += niterbicg

        if master
            println("Semi-implicit time step done")
        end

        #--------------------------------------------------------------
        #                           Output
        #--------------------------------------------------------------

        if output_steps
            if itime % noutput == 0
                iout = write_output(state, time, iout, cpu_start_time)
            end
        else
            if output
                iout = write_output(state, time, iout, cpu_start_time)
                output = false
                nextoutputtime = nextoutputtime + outputtimediff
                if nextoutputtime >= maxtime
                    nextoutputtime = maxtime
                end
            end
        end

        #-------------------------------------------
        #              Abort criteria
        #-------------------------------------------

        if !output_steps && time * tref >= maxtime
            if master
                naveragebicg = ntotalbicg / itime / 2

                println("")
                println(repeat("-", 80))
                println("Average Poisson iterations: ", naveragebicg)
                println(repeat("-", 80))
                println("")
            end

            break
        end
    end

    #-------------------------------------------
    #      Final output for output_steps
    #-------------------------------------------

    if output_steps
        if master
            naveragebicg = ntotalbicg / maxiter / 2

            println("")
            println(repeat("-", 80))
            println("Average Poisson iterations: ", naveragebicg)
            println(repeat("-", 80))
            println("")
        end
    end

    if master
        println("")
        println(repeat("-", 80))
        println(repeat(" ", 32), "PincFlow finished", repeat(" ", 33))
        println(repeat("-", 80))
    end

    return
end
