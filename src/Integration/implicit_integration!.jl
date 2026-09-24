"""
```julia
implicit_integration!(
    state::State,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    ntotalbicg::Integer,
    side::RHS,
    rayleigh_factor::AbstractFloat,
    iout::Integer,
    machine_start_time::DateTime,
)
```

Perform an implicit Euler step on the right-hand sides of the prognostic equations, solve the Poisson equation and correct the Exner-pressure, momentum and density fluctuations accordingly.

# Arguments

  - `state`: Model state.

  - `dtstage`: Fractional time step.

  - `time`: Simulation time.

  - `ntotalbicg`: BiCGSTAB-iterations counter.

  - `side`: Side of the equations.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

  - `iout`: Output counter.

  - `machine_start_time`: Wall-clock start time.
"""
function implicit_integration! end

function implicit_integration!(
    state::State,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    ntotalbicg::Integer,
    side::RHS,
    rayleigh_factor::AbstractFloat,
    iout::Integer,
    machine_start_time::DateTime,
)
    (; master) = state.domain

    modify_compressible_wind!(state, *)

    set_boundaries!(state, BoundaryPredictands())

    save_backups!(state, :w)

    update!(state, dtstage, U(), side, Implicit(), rayleigh_factor)
    update!(state, dtstage, V(), side, Implicit(), rayleigh_factor)
    update!(state, dtstage, W(), side, Implicit(), rayleigh_factor)

    set_boundaries!(state, BoundaryPredictands())

    update!(state, dtstage, RhoP(), side, Implicit(), rayleigh_factor)

    set_boundaries!(state, BoundaryPredictands())

    (errflagbicg, niterbicg) = apply_corrector!(state, dtstage, rayleigh_factor)

    if errflagbicg
        iout = write_output(state, time, iout, machine_start_time)
        error("BiCGSTAB errored! Output last state into record ", iout, ".")
    end

    ntotalbicg += niterbicg

    modify_compressible_wind!(state, /)

    set_boundaries!(state, BoundaryPredictands())

    return ntotalbicg
end
