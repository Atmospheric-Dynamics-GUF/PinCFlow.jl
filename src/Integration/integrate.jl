"""
    integrate(namelists::Namelists)

Main time integration loop for PinCFLow.

This function performs the complete time integration of the governing equations
using a semi-implicit time stepping scheme. It handles initialization,
time stepping, output, and finalization of the simulation.

# Arguments

  - `namelists::Namelists`: Configuration parameters for the simulation

# Algorithm Overview

The integration uses a 5-stage semi-implicit scheme:

 1. Explicit integration of LHS (advection/diffusion) over dt/2
 2. Implicit integration of RHS (pressure/acoustic terms) over dt/2
 3. Explicit integration of RHS over dt/2 (predictor step)
 4. Explicit integration of LHS over full dt
 5. Implicit integration of RHS over dt/2 (corrector step)

# Initialization

  - Sets up MPI communication and domain decomposition
  - Initializes all field variables and boundary conditions
  - Performs optional initial divergence cleaning
  - Initializes gravity wave ray tracing (MS-GWaM)
  - Handles restart from previous simulation if requested

# Time Stepping

Each time step involves:

  - CFL-based time step computation
  - Sponge layer application
  - Gravity wave propagation and saturation
  - Semi-implicit integration of momentum and continuity equations
  - Poisson solver for pressure correction
  - Boundary condition updates

# Output

  - Creates NetCDF output files
  - Writes field data at specified intervals
  - Supports both time-based and step-based output
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
