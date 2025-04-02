function apply_corrector!(
  state::State,
  dt::AbstractFloat,
  opt::AbstractIntegration,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; namelists, domain) = state
  (; model, zboundaries) = namelists.setting
  (; rhs) = state.poisson
  (; dpip) = state.variables.tendencies
  (; k0, k1) = domain

  # Initialize RHS and tolerance.
  rhs .= 0.0

  # Calculate RHS and tolreance reference.
  tolref = compute_rhs!(state, rhs, model)

  # Solve Poisson equation.
  (errflagbicg, niterbicg) =
    solve_poisson!(state, rhs, tolref, dt, opt, model, facray, facprs)

  # Return if an error occurred.
  if errflagbicg
    return (errflagbicg, niterbicg)
  end

  # Set boundaries of pressure correction.
  set_one_zonal_boundary_layer!(dpip, namelists, domain)
  set_one_meridional_boundary_layer!(dpip, namelists, domain)

  # Set vertical boundaries of dp.
  if zboundaries == SolidWallBoundaries()
    @views dpip[:, :, k0 - 1] .= dpip[:, :, k0]
    @views dpip[:, :, k1 + 1] .= dpip[:, :, k1]
  else
    error("Error in apply_corrector!: Unknown zboundaries!")
  end

  # Correct momentum and buoyancy.
  correct!(state, dt, opt, facray, facprs)

  # Return.
  return (errflagbicg, niterbicg)
end
