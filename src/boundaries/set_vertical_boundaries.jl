function set_vertical_boundaries!(
  state::State,
  variables::BoundaryPredictands,
  boundaries::SolidWallBoundaries,
)

  # Get all necessary fields.
  (; nbz) = state.namelists.domain
  (; k0, k1) = state.domain
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  # Set the density fluctuations at the boundaries to zero.
  for k in 1:nbz
    @views rho[:, :, k0 - k] .= -rho[:, :, k0 + k - 1]
    @views rho[:, :, k1 + k] .= -rho[:, :, k1 - k + 1]
  end

  # Set the density fluctuations at the boundaries to zero.
  for k in 1:nbz
    @views rhop[:, :, k0 - k] .= -rhop[:, :, k0 + k - 1]
    @views rhop[:, :, k1 + k] .= -rhop[:, :, k1 - k + 1]
  end

  # Set the vertical wind at the boundaries to zero.
  w[:, :, k0 - 1] .= 0.0
  w[:, :, k1] .= 0.0
  for k in 1:nbz
    @views w[:, :, k0 - k] .= -w[:, :, k0 + k - 2]
    @views w[:, :, k1 + k] .= -w[:, :, k1 - k]
  end

  # Set the horizontal-wind gradient at the boundaries to zero.
  for k in 1:nbz
    @views u[:, :, k0 - k] .= u[:, :, k0 + k - 1]
    @views u[:, :, k1 + k] .= u[:, :, k1 - k + 1]
    @views v[:, :, k0 - k] .= v[:, :, k0 + k - 1]
    @views v[:, :, k1 + k] .= v[:, :, k1 - k + 1]
  end

  # Set the pressure gradient at the boundaries to zero.
  @views pip[:, :, k0 - 1] .= pip[:, :, k0]
  @views pip[:, :, k1 + 1] .= pip[:, :, k1]

  return
end

function set_vertical_boundaries!(
  state::State,
  variables::BoundaryFluxes,
  boundaries::SolidWallBoundaries,
)

  # Get all necessary fields.
  (; k0, k1) = state.domain
  (; phirho, phirhop, phiu, phiv, phiw) = state.variables.fluxes

  # Set all vertical boundary fluxes to zero.

  # Density
  phirho[:, :, k0 - 1, 3] .= 0.0
  phirho[:, :, k1, 3] .= 0.0

  # Density fluctuations
  phirhop[:, :, k0 - 1, 3] .= 0.0
  phirhop[:, :, k1, 3] .= 0.0

  # Zonal momentum
  phiu[:, :, k0 - 1, 3] .= 0.0
  phiu[:, :, k1, 3] .= 0.0

  # Meridional momentum
  phiv[:, :, k0 - 1, 3] .= 0.0
  phiv[:, :, k1, 3] .= 0.0

  # Vertical momentum
  phiw[:, :, k0 - 2, 3] .= 0.0
  phiw[:, :, k1, 3] .= 0.0

  return
end

function set_vertical_boundaries!(
  state::State,
  variables::BoundaryGWIntegrals,
  boundaries::SolidWallBoundaries,
)
  # Get all necessary fields.
  (; k0, k1) = state.domain
  (; uu, uv, uw, vv, vw, etx, ety, utheta, vtheta, e) = state.wkb.integrals

  if steady_state || single_column
    @views uw[:, :, k0 - 1] .= uw[:, :, k0]
    @views uw[:, :, k1 + 1] .= uw[:, :, k1]
    @views vw[:, :, k0 - 1] .= vw[:, :, k0]
    @views vw[:, :, k1 + 1] .= vw[:, :, k1]
    @views e[:, :, k0 - 1] .= e[:, :, k0]
    @views e[:, :, k1 + 1] .= e[:, :, k1]
  else
    # no loop as Integrals also contains the tendencies
    @views uu[:, :, k0 - 1] .= uu[:, :, k0]
    @views uu[:, :, k1 + 1] .= uu[:, :, k1]
    @views uv[:, :, k0 - 1] .= uv[:, :, k0]
    @views uv[:, :, k1 + 1] .= uv[:, :, k1]
    @views uw[:, :, k0 - 1] .= uw[:, :, k0]
    @views uw[:, :, k1 + 1] .= uw[:, :, k1]
    @views vv[:, :, k0 - 1] .= vv[:, :, k0]
    @views vv[:, :, k1 + 1] .= vv[:, :, k1]
    @views vw[:, :, k0 - 1] .= vw[:, :, k0]
    @views vw[:, :, k1 + 1] .= vw[:, :, k1]
    @views etx[:, :, k0 - 1] .= etx[:, :, k0]
    @views etx[:, :, k1 + 1] .= etx[:, :, k1]
    @views ety[:, :, k0 - 1] .= ety[:, :, k0]
    @views ety[:, :, k1 + 1] .= ety[:, :, k1]
    @views utheta[:, :, k0 - 1] .= utheta[:, :, k0]
    @views utheta[:, :, k1 + 1] .= utheta[:, :, k1]
    @views vtheta[:, :, k0 - 1] .= vtheta[:, :, k0]
    @views vtheta[:, :, k1 + 1] .= vtheta[:, :, k1]
    @views e[:, :, k0 - 1] .= e[:, :, k0]
    @views e[:, :, k1 + 1] .= e[:, :, k1]
  end
  return
end

function set_vertical_boundaries!(
  state::State,
  variables::BoundaryGWTendencies,
  boundaries::SolidWallBoundaries,
)
  # Get all necessary fields.
  (; k0, k1) = state.domain
  (; dudt, dvdt, dthetadt) = state.wkb.integrals

  if steady_state || single_column
    @views dudt[:, :, k0 - 1] .= dudt[:, :, k0]
    @views dudt[:, :, k1 + 1] .= dudt[:, :, k1]
    @views dvdt[:, :, k0 - 1] .= dvdt[:, :, k0]
    @views dvdt[:, :, k1 + 1] .= dvdt[:, :, k1]
  else
    @views dudt[:, :, k0 - 1] .= dudt[:, :, k0]
    @views dudt[:, :, k1 + 1] .= dudt[:, :, k1]
    @views dvdt[:, :, k0 - 1] .= dvdt[:, :, k0]
    @views dvdt[:, :, k1 + 1] .= dvdt[:, :, k1]
    @views dthetadt[:, :, k0 - 1] .= dthetadt[:, :, k0]
    @views dthetadt[:, :, k1 + 1] .= dthetadt[:, :, k1]
  end
  return
end

function set_vertical_boundaries!(
  state::State,
  variables::BoundaryGWForces,
  boundaries::SolidWallBoundaries,
)
  # Get all necessary fields.
  (; k0, k1) = state.domain
  (; u, v, w) = state.wkb.gwmomforce

  @views u[:, :, k0 - 1] .= -u[:, :, k0]
  @views u[:, :, k1 + 1] .= -u[:, :, k1]
  @views v[:, :, k0 - 1] .= -v[:, :, k0]
  @views v[:, :, k1 + 1] .= -v[:, :, k1]
  @views w[:, :, k0 - 1] .= -w[:, :, k0]
  @views w[:, :, k1 + 1] .= -w[:, :, k1]
  return
end