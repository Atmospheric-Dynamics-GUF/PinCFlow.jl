function set_vertical_boundaries!(
  state::State,
  variables::BoundaryPredictands,
  boundaries::SolidWallBoundaries,
)

  # Get all necessary fields.
  (; nbz) = state.namelists.domain
  (; nz) = state.domain
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  # Set the density fluctuations at the boundaries to zero.
  for k in 1:(nbz)
    rho[:, :, -k + 1] = -rho[:, :, k]
    rho[:, :, nz + k] = -rho[:, :, nz - k + 1]
  end

  # Set the density fluctuations at the boundaries to zero.
  for k in 1:(nbz)
    rhop[:, :, -k + 1] = -rhop[:, :, k]
    rhop[:, :, nz + k] = -rhop[:, :, nz - k + 1]
  end

  # Set the vertical wind at the boundaries to zero.
  w[:, :, 0] .= 0.0
  w[:, :, nz] .= 0.0
  for k in 1:(nbz)
    w[:, :, -k] = -w[:, :, k]
    w[:, :, nz + k] = -w[:, :, nz - k]
  end

  # Set the horizontal-wind gradient at the boundaries to zero.
  for k in 1:(nbz)
    u[:, :, -k + 1] = u[:, :, k]
    u[:, :, nz + k] = u[:, :, nz - k + 1]
    v[:, :, -k + 1] = v[:, :, k]
    v[:, :, nz + k] = v[:, :, nz - k + 1]
  end

  # Set the pressure gradient at the boundaries to zero.
  pip[:, :, 0] = pip[:, :, 1]
  pip[:, :, nz + 1] = pip[:, :, nz]

  return
end

function set_vertical_boundaries!(
  state::State,
  variables::BoundaryFluxes,
  boundaries::SolidWallBoundaries,
)

  # Get all necessary fields.
  (; nz) = state.domain
  (; phirho, phirhop, phiu, phiv, phiw) = state.variables.fluxes

  # Set all vertical boundary fluxes to zero.

  # Density
  phirho[:, :, 0, 3] .= 0.0
  phirho[:, :, nz, 3] .= 0.0

  # Density fluctuations
  phirhop[:, :, 0, 3] .= 0.0
  phirhop[:, :, nz, 3] .= 0.0

  # Zonal momentum
  phiu[:, :, 0, 3] .= 0.0
  phiu[:, :, nz, 3] .= 0.0

  # Meridional momentum
  phiv[:, :, 0, 3] .= 0.0
  phiv[:, :, nz, 3] .= 0.0

  # Vertical momentum
  phiw[:, :, -1, 3] .= 0.0
  phiw[:, :, nz, 3] .= 0.0

  return
end
