"""
```julia
apply_corrector!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)::Tuple{Bool, <:Integer}
```

Perform the corrector step and return a tuple containing an error flag and the number of BicGStab iterations.

The right-hand side and the linear operator of the discrete Poisson equation are calculated. The equation is then solved for Exner-pressure differences, using a preconditioned BicGStab algorithm. Finally, the Exner-pressure, wind and density fluctuations are corrected accordingly.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

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
    rayleigh_factor::AbstractFloat,
)::Tuple{Bool, <:Integer}
    (; namelists, domain) = state
    (; model, zboundaries) = namelists.setting
    (; rhs) = state.poisson
    (; dpip) = state.variables.increments

    # Initialize RHS.
    rhs .= 0.0

    # Calculate RHS and tolreance reference.
    tolref = compute_rhs!(state, rhs, model)

    # Solve Poisson equation.
    (errflagbicg, niterbicg) =
        solve_poisson!(state, rhs, tolref, dt, rayleigh_factor)

    # Return if an error occurred.
    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    # Set boundaries of pressure correction.
    set_zonal_boundaries_of_field!(dpip, namelists, domain)
    set_meridional_boundaries_of_field!(dpip, namelists, domain)
    set_vertical_boundaries_of_field!(dpip, namelists, domain, zboundaries, +)

    # Correct momentum and buoyancy.
    correct!(state, dt, rayleigh_factor)

    return (errflagbicg, niterbicg)
end
