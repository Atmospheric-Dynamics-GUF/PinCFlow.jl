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

    # Save machine start time.
    machine_start_time = now()

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
    (; pstrattfc) = state.atmosphere
    (; p) = state.variables.predictands

    # Print information.
    if master
        println(repeat("-", 80))
        println(repeat(" ", 36), "PinCFlow")
        println(
            repeat(" ", 12),
            "developed by Rieper et al (2013) and Schmid et al (2021)",
        )
        println(repeat(" ", 28), "modified by many others")
        println(repeat("-", 80))
        println("")
        println("Date: ", Dates.Date(machine_start_time))
        println("Time: ", Dates.Time(machine_start_time))
        println("")
        println("Virtual topology: (nprocx, nprocy) = ", (nprocx, nprocy))
        println("")
    end

    #---------------------------------------------
    #        Initial divergence cleaning
    #---------------------------------------------

    if initialcleaning
        set_boundaries!(state, BoundaryPredictands())

        modify_compressible_wind!(state, *)
        set_boundaries!(state, BoundaryPredictands())

        (errflagbicg, niterbicg) = apply_corrector!(state, 1.0, 1.0, 1.0)

        modify_compressible_wind!(state, /)
        set_boundaries!(state, BoundaryPredictands())

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
            println("")
        end

        time = read_input!(state)

        if maxtime < time * tref
            error("Restart error: maxtime < time!")
        end

        set_boundaries!(state, BoundaryPredictands())

        if model == Compressible()
            pstrattfc .= p
            update_buoyancy_frequency!(state)
        end
    end

    #------------------------------------------
    #              Initial output
    #------------------------------------------

    # Create the output file.
    create_output(state)

    # Write the initial state.
    iout = write_output(state, time, iout, machine_start_time)

    output = false
    nextoutputtime = time * tref + outputtimediff

    #-----------------------------------------------------
    #                        Time loop
    #-----------------------------------------------------

    if master
        println("Starting the time loop...")
        println("")
    end

    if output_steps
        maxiterations = maxiter
    else
        maxiterations = 2^30
    end

    for itime in 1:maxiterations
        if master
            println(repeat("-", 80))
            println("Time step = ", itime)
            println("Time = ", time * tref, " seconds")
            println(repeat("-", 80))
            println("")
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
                    println("")
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
        #                           MS-GWaM
        #-----------------------------------------------------------------

        apply_saturation_scheme!(state, dt)

        for rkstage in 1:nstages
            propagate_rays!(state, dt, rkstage)
        end

        split_rays!(state)
        shift_rays!(state)
        merge_rays!(state)
        set_boundary_rays!(state)

        compute_mean_flow_effect!(state)

        #---------------------------------------------------------------
        #                   Semi-implicit time scheme
        #---------------------------------------------------------------

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
            println("")
            println("(1) Explicit integration of LHS over dt/2...")
            println("")
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

            update!(state, 0.5 * dt, rkstage, Rho())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                Rho(),
            )

            update!(state, 0.5 * dt, rkstage, RhoP(), LHS())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                RhoP(),
            )

            update!(state, 0.5 * dt, rkstage, P())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                P(),
            )

            # RK step for momentum

            set_boundaries!(state, BoundaryPredictands())

            update!(state, 0.5 * dt, rkstage, U(), LHS())
            update!(state, 0.5 * dt, rkstage, V(), LHS())
            update!(state, 0.5 * dt, rkstage, W(), LHS())
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                U(),
            )
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                V(),
            )
            apply_unified_sponge!(
                state,
                stepfrac[rkstage] * 0.5 * dt,
                time,
                W(),
            )
        end

        if model == Compressible()
            pstrattfc .= p
            apply_unified_sponge!(state, 0.5 * dt, time, PiP())
        end

        # (2) implicit integration of the linear right-hand sides of the
        #     equations for density fluctuations and momentum over half a
        #     time step, under consideration of the divergence constraint
        #     \psi^{n+1/2} = \psi^# + dt/2 Q(\psi^{n+1/2})

        if master
            println("(2) Implicit integration of RHS over dt/2...")
            println("")
        end

        modify_compressible_wind!(state, *)
        update_buoyancy_frequency!(state)

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
        (errflagbicg, niterbicg) = apply_corrector!(state, 0.5 * dt, 1.0, 1.0)

        if errflagbicg
            iout = write_output(state, time, iout, machine_start_time)
            if master
                println("Output last state into record", iout)
            end
            exit()
        end

        ntotalbicg += niterbicg

        modify_compressible_wind!(state, /)

        set_boundaries!(state, BoundaryPredictands())

        # put new state into var1 in order to save the advection velocities
        v1 = deepcopy(state.variables.predictands)

        # (3) explicit integration of the linear right-hand sides of the
        #     equations for density fluctuations and momentum over half a
        #     time step, under consideration of the divergence constraint
        #     \psi^\ast = \psi^n + dt/2 Q(\psi^n)

        if master
            println("(3) Explicit integration of RHS over dt/2...")
            println("")
        end

        # (3) uses updated pressure field and (5) adjusts pressure over half a
        # time step#
        state.variables.predictands.rho .= v0.rho
        state.variables.predictands.rhop .= v0.rhop
        state.variables.predictands.u .= v0.u
        state.variables.predictands.v .= v0.v
        state.variables.predictands.w .= v0.w

        if model == Compressible()
            state.variables.predictands.pip .= v0.pip
            state.variables.predictands.p .= v0.p
            update_buoyancy_frequency!(state)
            modify_compressible_wind!(state, *)
        end

        set_boundaries!(state, BoundaryPredictands())

        state.variables.backups.rhopold .= state.variables.predictands.rhop

        # update density fluctuations (rhopStar)
        update!(state, 0.5 * dt, RhoP(), RHS(), EXPL())

        # update winds (uStar, vStar, wStar)

        update!(state, 0.5 * dt, U(), RHS(), EXPL())
        update!(state, 0.5 * dt, V(), RHS(), EXPL())
        update!(state, 0.5 * dt, W(), RHS(), EXPL())

        # Update the Exner pressure.
        update!(state, 0.5 * dt, PiP())

        modify_compressible_wind!(state, /)

        set_boundaries!(state, BoundaryPredictands())

        # (4) explicit integration of convective and
        #     viscous-diffusive/turbulent fluxes over a full time step,
        #     with the advection velocity kept constant
        #     \psi^{\ast\ast} = \psi^\ast + A^dt (\psi^\ast, v^{n+1/2})

        if master
            println("(4) Explicit integration of LHS over dt...")
            println("")
        end

        v0 = deepcopy(v1)

        if model == Compressible()
            v1 = deepcopy(state.variables.predictands)
            pstrattfc .= v0.p
        end

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

            update!(state, dt, rkstage, Rho())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, Rho())

            update!(state, dt, rkstage, RhoP(), LHS())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, RhoP())

            update!(state, dt, rkstage, P())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, P())

            # RK step for momentum
            set_boundaries!(state, BoundaryPredictands())

            update!(state, dt, rkstage, U(), LHS())
            update!(state, dt, rkstage, V(), LHS())
            update!(state, dt, rkstage, W(), LHS())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, U())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, V())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, W())
        end

        if model == Compressible()
            pstrattfc .= p
            apply_unified_sponge!(state, dt, time, PiP())
        end

        # (5) implicit integration of the linear right-hand sides of the
        #     equations for density fluctuations and momentum over half a
        #     time step, under consideration of the divergence constraint
        #     \psi^{n+1} = \psi^{\ast\ast} + dt/2 Q(\psi^{n+1})

        if master
            println("(5) Implicit integration of RHS over dt/2...")
            println("")
        end

        modify_compressible_wind!(state, *)
        update_buoyancy_frequency!(state)

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
        (errflagbicg, niterbicg) = apply_corrector!(state, 0.5 * dt, 2.0, 1.0)

        modify_compressible_wind!(state, *)
        set_boundaries!(state, BoundaryPredictands())

        if errflagbicg
            iout = write_output(state, time, iout, machine_start_time)
            if master
                println("Output last state into record ", iout)
            end
            exit()
        end

        set_boundaries!(state, BoundaryPredictands())

        ntotalbicg += niterbicg

        if master
            println("...the semi-implicit time step is done.")
            println("")
        end

        #--------------------------------------------------------------
        #                           Output
        #--------------------------------------------------------------

        if output_steps
            if itime % noutput == 0
                iout = write_output(state, time, iout, machine_start_time)
            end
        else
            if output
                iout = write_output(state, time, iout, machine_start_time)
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

            println(repeat("-", 80))
            println("Average Poisson iterations: ", naveragebicg)
            println(repeat("-", 80))
            println("")
        end
    end

    if master
        println(repeat("-", 80))
        println(repeat(" ", 32), "PincFlow finished", repeat(" ", 33))
        println(repeat("-", 80))
    end

    return
end
