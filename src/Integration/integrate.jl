"""
```julia
integrate(namelists::Namelists)
```

Initialize the model state and integrate it in time.

This method performs the complete time integration of the governing equations,
using a semi-implicit time stepping scheme. It handles initialization,
time stepping and output of the simulation data.

The initialization begins with the construction of the model state (an instance of the composite type `State`), which involves the setup of the MPI parallelization and the definition of all arrays that are needed repeatedly during the simulation. This is followed by an (optional) initial cleaning, in which the Poisson solver is called to ensure that the initial dynamic fields satisfy the divergence constraint imposed by the thermodynamic energy equation. Afterwards, the initialization of MS-GWaM is completed by adding ray volumes to the previously defined arrays. If the simulation is supposed to start from a previous model state, the fields are then overwritten with the data in the corresponding input file. Finally, the output file is created and the initial state is written into it.

At the beginning of each time-loop iteration, the time step is determined from several stability criteria, using `compute_time_step`. In case the updated simulation time is later than the next output time, the time step is corrected accordingly. Subsequently, the damping coefficients of the sponges (which may depend on the time step) are calculated. Following this, MS-GWaM updates the unresolved gravity-wave field and computes the corresponding mean-flow impact. Afterwards, the resolved flow is updated in a semi-implicit time step, comprised of the following stages.

 1. Explicit RK3 integration of LHS over ``\\Delta t / 2``.

 1. Implicit Euler integration of RHS over ``\\Delta t / 2``.

 1. Explicit Euler integration of RHS over ``\\Delta t / 2``.

 1. Explicit RK3 integration of LHS over ``\\Delta t``.

 1. Implicit Euler integration of RHS over ``\\Delta t / 2``.

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

  - [`PinCFlow.Update.compute_sponges!`](@ref)

  - [`PinCFlow.Integration.wkb_integration!`](@ref)

  - [`PinCFlow.Integration.synchronize_density_fluctuations!`](@ref)

  - [`PinCFlow.Integration.explicit_integration!`](@ref)

  - [`PinCFlow.Integration.implicit_integration!`](@ref)

  - [`PinCFlow.Integration.reset_predictands!`](@ref)
"""
function integrate end

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

    (; npx, npy, npz) = state.namelists.domain
    (; initial_cleaning) = state.namelists.poisson
    (; dtmin) = state.namelists.discretization
    (; restart, tmax, output_interval, output_steps, iterations, nout) =
        state.namelists.output
    (; tref) = state.constants
    (; master) = state.domain

    # Print information.
    if master
        println(repeat("-", 80))
        println(repeat(" ", 34), "PinCFlow.jl")
        println("")
        println(repeat(" ", 34), "developed by")
        println(repeat(" ", 30), "Rieper et al. (2013)")
        println(repeat(" ", 28), "Muraschko et al. (2014)")
        println(repeat(" ", 29), "Boeloeni et al. (2016)")
        println(repeat(" ", 29), "Wilhelm et al. (2018)")
        println(repeat(" ", 31), "Wei et al. (2019)")
        println(repeat(" ", 30), "Schmid et al. (2021)")
        println(repeat(" ", 30), "Jochum et al. (2025)")
        println("")
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

    if initial_cleaning
        modify_compressible_wind!(state, *)

        set_boundaries!(state, BoundaryPredictands())

        (errflagbicg, niterbicg) = apply_corrector!(state, 1.0, 1.0)

        if errflagbicg
            create_output(state, machine_start_time)
            iout = write_output(state, time, iout, machine_start_time)
            if master
                println("Output last state into record ", iout, ".")
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

        if tmax < time * tref
            error("Restart error: tmax < time!")
        end

        set_boundaries!(state, BoundaryPredictands())

        synchronize_compressible_atmosphere!(state, state.variables.predictands)
    end

    #------------------------------------------
    #              Initial output
    #------------------------------------------

    # Create the output file.
    create_output(state, machine_start_time)

    # Write the initial state.
    iout = write_output(state, time, iout, machine_start_time)

    # Prepare the next output.
    output = false
    nextoutputtime = time * tref + output_interval

    #-----------------------------------------------------
    #                        Time loop
    #-----------------------------------------------------

    if master
        println("Starting the time loop...")
        println("")
    end

    if output_steps
        maxiterations = iterations
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
            if (time + dt) * tref + dtmin > nextoutputtime
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
        #                           Sponges
        #-----------------------------------------------------------------

        compute_sponges!(state, dt, time)

        #-----------------------------------------------------------------
        #                           MS-GWaM
        #-----------------------------------------------------------------

        wkb_integration!(state, dt)

        #-----------------------------------------------------------------
        #                         Turbulence 
        #-----------------------------------------------------------------

        compute_turbulence_diffusion!(state)

        save_backups!(state, :u, :v)

        turbulent_diffusion!(state, dt)

        set_boundaries!(state, BoundaryPredictands())

        #---------------------------------------------------------------
        #                   Semi-implicit time scheme
        #---------------------------------------------------------------

        if master
            println("Beginning a semi-implicit time step...")
            println("")
        end

        synchronize_density_fluctuations!(state)

        set_boundaries!(state, BoundaryPredictands())

        (p0, chi0) = backup_predictands(state)

        compute_fluxes!(state, p0, Theta())

        if master
            println("(1) Explicit integration of LHS over dt/2...")
            println("")
        end

        explicit_integration!(state, p0, 0.5 * dt, time, LHS())

        if master
            println("(2) Implicit integration of RHS over dt/2...")
            println("")
        end

        ntotalbicg = implicit_integration!(
            state,
            0.5 * dt,
            time,
            ntotalbicg,
            RHS(),
            1.0,
            iout,
            machine_start_time,
        )

        p1 = deepcopy(state.variables.predictands)

        if master
            println("(3) Explicit integration of RHS over dt/2...")
            println("")
        end

        reset_predictands!(state, p0, chi0)

        explicit_integration!(state, p0, 0.5 * dt, time, RHS())

        if master
            println("(4) Explicit integration of LHS over dt...")
            println("")
        end

        p0 = deepcopy(p1)

        synchronize_compressible_atmosphere!(state, p0)

        turbulence_integration!(state, p0, dt)

        explicit_integration!(state, p0, dt, time, LHS())

        if master
            println("(5) Implicit integration of RHS over dt/2...")
            println("")
        end

        ntotalbicg = implicit_integration!(
            state,
            0.5 * dt,
            time,
            ntotalbicg,
            RHS(),
            2.0,
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
            if itime % nout == 0
                iout = write_output(state, time, iout, machine_start_time)
            end
        else
            if output
                iout = write_output(state, time, iout, machine_start_time)
                output = false
                nextoutputtime = nextoutputtime + output_interval
                if nextoutputtime >= tmax
                    nextoutputtime = tmax
                end
            end
        end

        #-------------------------------------------
        #              Abort criteria
        #-------------------------------------------

        if !output_steps && time * tref >= tmax
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
            naveragebicg = ntotalbicg / iterations / 2

            println(repeat("-", 80))
            println("Average Poisson iterations: ", naveragebicg)
            println(repeat("-", 80))
            println("")
        end
    end

    if master
        println(repeat("-", 80))
        println(repeat(" ", 30), "PinCFlow.jl finished")
        println(repeat("-", 80))
    end

    return
end
