function apply_corrector!(
  state::State,
  dt::AbstractFloat,
  errflagbicg::Bool,
  niter::Integer,
  opt::AbstractIntegration,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; namelists, domain) = state
  (; model) = namelists.setting
  (; zboundaries) = namelists.boundaries
  (; nz) = domain
  (; rhs) = state.poisson
  (; dpip) = state.variables.tendencies

  # Initialize RHS and tolerance.
  rhs .= 0.0
  tolref = 0.0

  # Calculate RHS.
  compute_rhs!(state, rhs, tolref, dt, model)

  # Solve Poisson equation.
  solve_poisson!(
    state,
    rhs,
    tolref,
    dt,
    errflagbicg,
    niter,
    opt,
    model,
    facray,
    facprs,
  )

  # Return if an error occurred.
  if errflagbicg
    return
  end

  # Set boundaries of pressure correction.
  set_zonal_boundaries_of_field!(dpip, namelists, domain)
  set_meridional_boundaries_of_field!(dpip, namelists, domain)

  # Set vertical boundaries of dp.
  if zboundaries == SolidWallBoundaries()
    dpip[:, :, 0] = dpip[:, :, 1]
    dpip[:, :, nz + 1] = dpip[:, :, nz]
  else
    error("Error in apply_corrector: Unknown zboundaries!")
  end

  # Correct momentum and buoyancy.
  correct!(state, dt, opt, facray, facprs)

  # Return.
  return
end
