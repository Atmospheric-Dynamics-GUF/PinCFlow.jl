"""
```julia
integrate(namelists::Namelists)
```

Initialize the model state and integrate it in time.

This method performs the complete time integration of the governing equations,
using a semi-implicit time stepping scheme. It handles initialization,
time stepping and output of the simulation data.

The initialization begins with the construction of the model state (an instance of the composite type `State`), which involves the setup of the MPI parallelization and the definition of all arrays that are needed repeatedly during the simulation. This is followed by an (optional) initial cleaning, in which the Poisson solver is called to ensure that the initial dynamic fields satisfy the divergence constraint imposed by the thermodynamic energy equation. Afterwards, the initialization of MS-GWaM is completed by adding ray volumes to the previously defined arrays. If the simulation is supposed to start from a previous model state, the fields are then overwritten with the data in the corresponding input file. Finally, the output file is created and the initial state is written into it.

At the beginning of each time-loop iteration, the time step is determined from several stability criteria, using `compute_time_step`. In case the updated simulation time is later than the next output time, the time step is corrected accordingly. Subsequently, the damping coefficient of the sponge layer (which may depend on the time step) is calculated. Following this, MS-GWaM updates the unresolved gravity-wave field and computes the corresponding mean-flow impact. Afterwards, the resolved flow is updated in a semi-implicit time step, comprised of the following stages.

 1. Explicit RK3 integration of LHS over ``\\Delta t / 2``.
 2. Implicit Euler integration of RHS over ``\\Delta t / 2``.
 3. Explicit Euler integration of RHS over ``\\Delta t / 2``.
 4. Explicit RK3 integration of LHS over ``\\Delta t``.
 5. Implicit Euler integration of RHS over ``\\Delta t / 2``.

Therein, the left-hand sides of the equations include advective fluxes, diffusion terms, rotation and heating, whereas the pressure gradient, buoyancy term and momentum-flux divergence due to unresolved gravity waves are on the right-hand sides. Boundary conditions are enforced continuously. At the end of the time step, the updated fields are written into the output file if the next output time has been reached.

# Arguments

  - `namelists`: Namelists with all model parameters.

# See also

  - [`PinCFlow.Types.State`](@ref)
  - [`PinCFlow.Integration.modify_compressible_wind!`](@ref)
  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)
  - [`PinCFlow.PoissonSolver.apply_corrector!`](@ref)
  - [`PinCFlow.Output.create_output`](@ref)
  - [`PinCFlow.Output.write_output`](@ref)
  - [`PinCFlow.MSGWaM.RayUpdate.initialize_rays!`](@ref)
  - [`PinCFlow.Output.read_input!`](@ref)
  - [`PinCFlow.Integration.synchronize_compressible_atmosphere!`](@ref)
  - [`PinCFlow.Integration.compute_time_step`](@ref)
  - [`PinCFlow.Update.compute_sponge!`](@ref)
  - [`PinCFlow.MSGWaM.RayUpdate.apply_saturation_scheme!`](@ref)
  - [`PinCFlow.MSGWaM.RayUpdate.propagate_rays!`](@ref)
  - [`PinCFlow.MSGWaM.RayUpdate.split_rays!`](@ref)
  - [`PinCFlow.MSGWaM.RayUpdate.shift_rays!`](@ref)
  - [`PinCFlow.MSGWaM.RayUpdate.merge_rays!`](@ref)
  - [`PinCFlow.MSGWaM.BoundaryRays.set_boundary_rays!`](@ref)
  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_mean_flow_effect!`](@ref)
  - [`PinCFlow.Integration.synchronize_density_fluctuations!`](@ref)
  - [`PinCFlow.FluxCalculator.reconstruct!`](@ref)
  - [`PinCFlow.FluxCalculator.compute_fluxes!`](@ref)
  - [`PinCFlow.Update.update!`](@ref)
  - [`PinCFlow.Update.apply_unified_sponge!`](@ref)
  - [`PinCFlow.Integration.save_backups!`](@ref)
  - [`PinCFlow.Integration.reset_predictands!`](@ref)
"""
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
    (; npx, npy, npz) = state.namelists.domain
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
        println("Virtual topology: (npx, npy, npz) = ", (npx, npy, npz))
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
            create_output(state)
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

        wkb_integration!(state, dt)

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

        explicit_integration!(state, p0, 0.5 * dt, time, LHS())

        if master
            println("(2) Implicit integration of RHS over dt/2...")
            println("")
        end

        implicit_integration!(
            state,
            0.5 * dt,
            time,
            ntotalbicg,
            RHS(),
            iout,
            machine_start_time,
        )

        p1 = deepcopy(state.variables.predictands)

        if master
            println("(3) Explicit integration of RHS over dt/2...")
            println("")
        end

        reset_predictands!(state, p0)

        explicit_integration!(state, p0, 0.5 * dt, time, RHS())

        if master
            println("(4) Explicit integration of LHS over dt...")
            println("")
        end

        p0 = deepcopy(p1)

        synchronize_compressible_atmosphere!(state, p0)

        explicit_integration!(state, p0, dt, time, LHS())

        if master
            println("(5) Implicit integration of RHS over dt/2...")
            println("")
        end

        implicit_integration!(
            state,
            0.5 * dt,
            time,
            ntotalbicg,
            RHS(),
            iout,
            machine_start_time,
        )

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
