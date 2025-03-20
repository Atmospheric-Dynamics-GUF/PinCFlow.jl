function reconstruct!(state::State)
  reconstruct!(state, Rho())
  reconstruct!(state, RhoP())
  reconstruct!(state, U())
  reconstruct!(state, V())
  reconstruct!(state, W())
  return
end

function reconstruct!(state::State, variable::Rho)
  (; nbx, nby, nbz) = state.namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = state.domain
  (; rho) = state.variables.predictands
  (; rhobar) = state.variables.auxiliaries
  (; rhotilde) = state.variables.reconstructions
  (; pstrattfc) = state.atmosphere

  for ix in (-nbx):(nx + nbx)
    for jy in (-nby):(ny + nby)
      for kz in 0:(nz + 1)
        rhobar[ix, jy, kz] = rho[ix, jy, kz] / pstrattfc[ix, jy, kz]
      end
    end
  end
  apply_3d_muscl!(
    rhobar,
    rhotilde,
    nx + 2 * nbx + 1,
    ny + 2 * nby + 1,
    nz + 2 * nbz + 1,
    limitertype,
  )

  return
end

function reconstruct!(state::State, variable::RhoP)
  (; nbx, nby, nbz) = state.namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = state.domain
  (; rhop) = state.variables.predictands
  (; rhopbar) = state.variables.auxiliaries
  (; rhoptilde) = state.variables.reconstructions
  (; pstrattfc) = state.atmosphere

  for ix in (-nbx):(nx + nbx)
    for jy in (-nby):(ny + nby)
      for kz in 0:(nz + 1)
        rhopbar[ix, jy, kz] = rhop[ix, jy, kz] / pstrattfc[ix, jy, kz]
      end
    end
  end
  apply_3d_muscl!(
    rhopbar,
    rhoptilde,
    nx + 2 * nbx + 1,
    ny + 2 * nby + 1,
    nz + 2 * nbz + 1,
    limitertype,
  )

  return
end

function reconstruct!(state::State, variable::U)
  (; nbx, nby, nbz) = state.namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = state.domain
  (; rho, u) = state.variables.predictands
  (; ubar) = state.variables.auxiliaries
  (; utilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  for ix in (-nbx):(nx + nbx - 1)
    for jy in (-nby):(ny + nby)
      for kz in 0:(nz + 1)
        rhoedge =
          0.5 * (
            rho[ix, jy, kz] +
            rho[ix + 1, jy, kz] +
            rhostrattfc[ix, jy, kz] +
            rhostrattfc[ix + 1, jy, kz]
          )
        pedge = 0.5 * (pstrattfc[ix, jy, kz] + pstrattfc[ix + 1, jy, kz])
        ubar[ix, jy, kz] = u[ix, jy, kz] * rhoedge / pedge
      end
    end
  end

  apply_3d_muscl!(
    ubar,
    utilde,
    nx + 2 * nbx + 1,
    ny + 2 * nby + 1,
    nz + 2 * nbz + 1,
    limitertype,
  )

  return
end

function reconstruct!(state::State, variable::V)
  (; nbx, nby, nbz) = state.namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = state.domain
  (; rho, v) = state.variables.predictands
  (; vbar) = state.variables.auxiliaries
  (; vtilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  for ix in (-nbx):(nx + nbx)
    for jy in (-nby):(ny + nby - 1)
      for kz in 0:(nz + 1)
        rhoedge =
          0.5 * (
            rho[ix, jy, kz] +
            rho[ix, jy + 1, kz] +
            rhostrattfc[ix, jy, kz] +
            rhostrattfc[ix, jy + 1, kz]
          )
        pedge = 0.5 * (pstrattfc[ix, jy, kz] + pstrattfc[ix, jy + 1, kz])
        vbar[ix, jy, kz] = v[ix, jy, kz] * rhoedge / pedge
      end
    end
  end

  apply_3d_muscl!(
    vbar,
    vtilde,
    nx + 2 * nbx + 1,
    ny + 2 * nby + 1,
    nz + 2 * nbz + 1,
    limitertype,
  )

  return
end

function reconstruct!(state::State, variable::W)
  (; namelists, domain, grid) = state
  (; nbx, nby, nbz, nprocx, nprocy) = namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = domain
  (; jac) = grid
  (; predictands) = state.variables
  (; rho, w) = predictands
  (; wbar) = state.variables.auxiliaries
  (; wtilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  wbar[:, :, 0:(nz + 1)] = w[:, :, 0:(nz + 1)]
  for ix in 1:nx
    for jy in 1:ny
      for kz in 0:(nz + 1)
        wbar[ix, jy, kz] = compute_vertical_wind(ix, jy, kz, predictands, grid)
      end
    end
  end
  set_zonal_boundaries_of_field!(wbar, namelists, domain)
  set_meridional_boundaries_of_field!(wbar, namelists, domain)
  for ix in (-nbx):(nx + nbx)
    for jy in (-nby):(ny + nby)
      for kz in 0:(nz + 1)
        rhoedgeu =
          (
            jac[ix, jy, kz + 1] * (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
            jac[ix, jy, kz] *
            (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
          ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
        pedgeu =
          (
            jac[ix, jy, kz + 1] * pstrattfc[ix, jy, kz] +
            jac[ix, jy, kz] * pstrattfc[ix, jy, kz + 1]
          ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
        wbar[ix, jy, kz] *= rhoedgeu / pedgeu
      end
    end
  end

  apply_3d_muscl!(
    wbar,
    wtilde,
    nx + 2 * nbx + 1,
    ny + 2 * nby + 1,
    nz + 2 * nbz + 1,
    limitertype,
  )

  return
end
