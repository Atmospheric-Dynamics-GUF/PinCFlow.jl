function apply_corrector!(
    state::State,
    dt::AbstractFloat,
    facray::AbstractFloat,
    facprs::AbstractFloat,
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
    (errflagbicg, niterbicg) =
        solve_poisson!(state, rhs, tolref, dt, facray, facprs)

    # Return if an error occurred.
    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    # Set boundaries of pressure correction.
    set_zonal_boundaries_of_field!(dpip, namelists, domain)
    set_meridional_boundaries_of_field!(dpip, namelists, domain)
    set_vertical_boundaries_of_field!(dpip, namelists, domain, zboundaries, +)

    # Correct momentum and buoyancy.
    correct!(state, dt, facray, facprs)

    # Return.
    return (errflagbicg, niterbicg)
end
