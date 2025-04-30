function integrate(namelists::Namelists)

    #-------------------------------------------------
    #                     Setup
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
    (; initialcleaning) = state.namelists.poisson
    (; dtmin_dim) = state.namelists.discretization
    (; restart, maxtime, outputtimediff, output_steps, maxiter, noutput) =
        state.namelists.output
    (; nstages, stepfrac) = state.time
    (; tref) = state.constants
    (; master) = state.domain

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
        modify_compressible_wind!(state, *)

        set_boundaries!(state, BoundaryPredictands())

        (errflagbicg, niterbicg) = apply_corrector!(state, 1.0, 1.0, 1.0)

        if errflagbicg
            iout = write_output(state, time, iout, machine_start_time)
            if master
                println("Output last state into record", iout)
            end
            exit()
        end

        modify_compressible_wind!(state, /)

        set_boundaries!(state, BoundaryPredictands())
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

        synchronize_compressible_atmosphere!(state, state.variables.predictands)
    end

    #------------------------------------------
    #              Initial output
    #------------------------------------------

    # Create the output file.
    create_output(state)

    # Write the initial state.
    iout = write_output(state, time, iout, machine_start_time)

    # Prepare the next output.
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

        # Correct dt to hit desired output time.
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

        if master
            println("Beginning a semi-implicit time step...")
            println("")
        end

        synchronize_density_fluctuations!(state)

        set_boundaries!(state, BoundaryPredictands())

        p0 = deepcopy(state.variables.predictands)

        if master
            println("(1) Explicit integration of LHS over dt/2...")
            println("")
        end

        for rkstage in 1:nstages
            reconstruct!(state)
            set_boundaries!(state, BoundaryReconstructions())

            compute_fluxes!(state, p0)
            set_boundaries!(state, BoundaryFluxes())

            save_backups!(state, :rho)

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

            set_boundaries!(state, BoundaryPredictands())
        end

        synchronize_compressible_atmosphere!(state, state.variables.predictands)

        apply_unified_sponge!(state, 0.5 * dt, time, PiP())

        if master
            println("(2) Implicit integration of RHS over dt/2...")
            println("")
        end

        modify_compressible_wind!(state, *)

        set_boundaries!(state, BoundaryPredictands())

        save_backups!(state, :w)

        update!(state, 0.5 * dt, U(), RHS(), Implicit(), 1.0)
        update!(state, 0.5 * dt, V(), RHS(), Implicit(), 1.0)
        update!(state, 0.5 * dt, W(), RHS(), Implicit(), 1.0)

        set_boundaries!(state, BoundaryPredictands())

        update!(state, 0.5 * dt, RhoP(), RHS(), Implicit(), 1.0)

        set_boundaries!(state, BoundaryPredictands())

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

        p1 = deepcopy(state.variables.predictands)

        if master
            println("(3) Explicit integration of RHS over dt/2...")
            println("")
        end

        reset_predictands!(state, p0)

        synchronize_compressible_atmosphere!(state, state.variables.predictands)

        modify_compressible_wind!(state, *)

        set_boundaries!(state, BoundaryPredictands())

        save_backups!(state, :rhop, :u, :v, :w)

        update!(state, 0.5 * dt, RhoP(), RHS(), Explicit())

        update!(state, 0.5 * dt, U(), RHS(), Explicit())
        update!(state, 0.5 * dt, V(), RHS(), Explicit())
        update!(state, 0.5 * dt, W(), RHS(), Explicit())

        update!(state, 0.5 * dt, PiP())

        modify_compressible_wind!(state, /)

        set_boundaries!(state, BoundaryPredictands())

        if master
            println("(4) Explicit integration of LHS over dt...")
            println("")
        end

        p0 = deepcopy(p1)

        # I'm not sure why we use P at dt/2 here, instead of the initial P...
        synchronize_compressible_atmosphere!(state, p0)

        for rkstage in 1:nstages
            reconstruct!(state)
            set_boundaries!(state, BoundaryReconstructions())

            compute_fluxes!(state, p0)
            set_boundaries!(state, BoundaryFluxes())

            save_backups!(state, :rho)

            update!(state, dt, rkstage, Rho())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, Rho())

            update!(state, dt, rkstage, RhoP(), LHS())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, RhoP())

            update!(state, dt, rkstage, P())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, P())

            set_boundaries!(state, BoundaryPredictands())

            update!(state, dt, rkstage, U(), LHS())
            update!(state, dt, rkstage, V(), LHS())
            update!(state, dt, rkstage, W(), LHS())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, U())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, V())
            apply_unified_sponge!(state, stepfrac[rkstage] * dt, time, W())

            set_boundaries!(state, BoundaryPredictands())
        end

        synchronize_compressible_atmosphere!(state, state.variables.predictands)

        apply_unified_sponge!(state, dt, time, PiP())

        if master
            println("(5) Implicit integration of RHS over dt/2...")
            println("")
        end

        modify_compressible_wind!(state, *)

        set_boundaries!(state, BoundaryPredictands())

        save_backups!(state, :w)

        update!(state, 0.5 * dt, U(), RHS(), Implicit(), 2.0)
        update!(state, 0.5 * dt, V(), RHS(), Implicit(), 2.0)
        update!(state, 0.5 * dt, W(), RHS(), Implicit(), 2.0)

        set_boundaries!(state, BoundaryPredictands())

        update!(state, 0.5 * dt, RhoP(), RHS(), Implicit(), 2.0)

        set_boundaries!(state, BoundaryPredictands())

        (errflagbicg, niterbicg) = apply_corrector!(state, 0.5 * dt, 2.0, 1.0)

        if errflagbicg
            iout = write_output(state, time, iout, machine_start_time)
            if master
                println("Output last state into record ", iout)
            end
            exit()
        end

        ntotalbicg += niterbicg

        modify_compressible_wind!(state, /)

        set_boundaries!(state, BoundaryPredictands())

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
