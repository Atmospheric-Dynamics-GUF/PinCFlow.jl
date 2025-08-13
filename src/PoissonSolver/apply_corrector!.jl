"""
```julia
apply_corrector!(
    state::State,
    dt::AbstractFloat,
    facray::AbstractFloat,
)
```

Perform the corrector step by computing the right-hand side and linear operator of the discrete Poisson equation, solving it and correcting the Exner-pressure, wind and buoyancy accordingly.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `facray`: Factor by which the Rayleigh-damping coefficient is multiplied.

# Returns

  - `::Bool`: Error flag.

  - `::Integer`: Number of iterations.

# See also

  - [`PinCFlow.PoissonSolver.compute_rhs!`](@ref)

  - [`PinCFlow.PoissonSolver.solve_poisson!`](@ref)

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)

  - [`PinCFlow.PoissonSolver.correct!`](@ref)
"""
function apply_corrector! end

function apply_corrector!(
    state::State,
    dt::AbstractFloat,
    facray::AbstractFloat,
)
    (; namelists, domain) = state
    (; model, zboundaries) = namelists.setting
    (; rhs) = state.poisson
    (; dpip) = state.variables.tendencies

    # Initialize RHS.
    rhs .= 0.0

    # Calculate RHS and tolreance reference.
    tolref = compute_rhs!(state, rhs, model)

    # Solve Poisson equation.
    (errflagbicg, niterbicg) = solve_poisson!(state, rhs, tolref, dt, facray)

    # Return if an error occurred.
    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    # Set boundaries of pressure correction.
    set_zonal_boundaries_of_field!(dpip, namelists, domain)
    set_meridional_boundaries_of_field!(dpip, namelists, domain)
    set_vertical_boundaries_of_field!(dpip, namelists, domain, zboundaries, +)

    # Correct momentum and buoyancy.
    correct!(state, dt, facray)

    return (errflagbicg, niterbicg)
end
