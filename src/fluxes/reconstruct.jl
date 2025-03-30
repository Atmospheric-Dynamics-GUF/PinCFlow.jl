function reconstruct!(state::State)
  reconstruct!(state, Rho())
  reconstruct!(state, RhoP())
  reconstruct!(state, U())
  reconstruct!(state, V())
  reconstruct!(state, W())
  return
end

function reconstruct!(state::State, variable::Rho)
  (; limitertype) = state.namelists.discretization
  (; k0, k1, nxx, nyy, nzz) = state.domain
  (; rho) = state.variables.predictands
  (; phi) = state.variables.auxiliaries
  (; rhotilde) = state.variables.reconstructions
  (; pstrattfc) = state.atmosphere

  for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
    phi[ix, jy, kz] = rho[ix, jy, kz] / pstrattfc[ix, jy, kz]
  end
  apply_3d_muscl!(phi, rhotilde, nxx, nyy, nzz, limitertype)

  return
end

function reconstruct!(state::State, variable::RhoP)
  (; limitertype) = state.namelists.discretization
  (; k0, k1, nxx, nyy, nzz) = state.domain
  (; rhop) = state.variables.predictands
  (; phi) = state.variables.auxiliaries
  (; rhoptilde) = state.variables.reconstructions
  (; pstrattfc) = state.atmosphere

  for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
    phi[ix, jy, kz] = rhop[ix, jy, kz] / pstrattfc[ix, jy, kz]
  end
  apply_3d_muscl!(phi, rhoptilde, nxx, nyy, nzz, limitertype)

  return
end

function reconstruct!(state::State, variable::U)
  (; limitertype) = state.namelists.discretization
  (; k0, k1, nxx, nyy, nzz) = state.domain
  (; rho, u) = state.variables.predictands
  (; phi) = state.variables.auxiliaries
  (; utilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:(nxx - 1)
    rhoedge =
      0.5 * (
        rho[ix, jy, kz] +
        rho[ix + 1, jy, kz] +
        rhostrattfc[ix, jy, kz] +
        rhostrattfc[ix + 1, jy, kz]
      )
    pedge = 0.5 * (pstrattfc[ix, jy, kz] + pstrattfc[ix + 1, jy, kz])
    phi[ix, jy, kz] = u[ix, jy, kz] * rhoedge / pedge
  end

  apply_3d_muscl!(phi, utilde, nxx, nyy, nzz, limitertype)

  return
end

function reconstruct!(state::State, variable::V)
  (; limitertype) = state.namelists.discretization
  (; k0, k1, nxx, nyy, nzz) = state.domain
  (; rho, v) = state.variables.predictands
  (; phi) = state.variables.auxiliaries
  (; vtilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  for kz in (k0 - 1):(k1 + 1), jy in 1:(nyy - 1), ix in 1:nxx
    rhoedge =
      0.5 * (
        rho[ix, jy, kz] +
        rho[ix, jy + 1, kz] +
        rhostrattfc[ix, jy, kz] +
        rhostrattfc[ix, jy + 1, kz]
      )
    pedge = 0.5 * (pstrattfc[ix, jy, kz] + pstrattfc[ix, jy + 1, kz])
    phi[ix, jy, kz] = v[ix, jy, kz] * rhoedge / pedge
  end

  apply_3d_muscl!(phi, vtilde, nxx, nyy, nzz, limitertype)

  return
end

function reconstruct!(state::State, variable::W)
  (; namelists, domain, grid) = state
  (; limitertype) = state.namelists.discretization
  (; i0, i1, j0, j1, k0, k1, nxx, nyy, nzz) = domain
  (; jac) = grid
  (; predictands) = state.variables
  (; rho, w) = predictands
  (; phi) = state.variables.auxiliaries
  (; wtilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  @views phi[:, :, (k0 - 1):(k1 + 1)] .= w[:, :, (k0 - 1):(k1 + 1)]
  for kz in (k0 - 1):(k1 + 1), jy in j0:j1, ix in i0:i1
    phi[ix, jy, kz] = compute_vertical_wind(ix, jy, kz, predictands, grid)
  end
  set_zonal_boundaries_of_field!(phi, namelists, domain)
  set_meridional_boundaries_of_field!(phi, namelists, domain)
  for kz in (k0 - 1):(k1 + 1), jy in 1:nyy, ix in 1:nxx
    rhoedgeu =
      (
        jac[ix, jy, kz + 1] * (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
        jac[ix, jy, kz] * (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
      ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
    pedgeu =
      (
        jac[ix, jy, kz + 1] * pstrattfc[ix, jy, kz] +
        jac[ix, jy, kz] * pstrattfc[ix, jy, kz + 1]
      ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
    phi[ix, jy, kz] *= rhoedgeu / pedgeu
  end

  apply_3d_muscl!(phi, wtilde, nxx, nyy, nzz, limitertype)

  return
end
