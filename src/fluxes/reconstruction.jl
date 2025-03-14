function apply_1d_muscl!(
  phi::Vector{<:AbstractFloat},
  phitilde::Matrix{<:AbstractFloat},
  phisize::Integer,
)

  # Initialize phitilde.
  phitilde .= 1000.0

  for i in 2:(phisize - 1)
    deltal = phi[i] - phi[i - 1]
    deltar = phi[i + 1] - phi[i]

    if deltar == 0.0
      phitilde[i, 2] = phi[i]
      phitilde[i, 1] = phi[i]
    else
      if deltal == 0.0
        theta = deltal / deltar
        s = (2.0 + theta) / 3.0
        sigmal = max(0.0, min(2 * theta, s, 2.0))

        phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
        phitilde[i, 1] = phi[i]
      else
        theta = deltal / deltar

        s = (2.0 + theta) / 3.0
        sigmal = max(0.0, min(2 * theta, s, 2.0))

        s = (2.0 + 1.0 / theta) / 3.0
        sigmar = max(0.0, min(2 / theta, s, 2.0))

        phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
        phitilde[i, 1] = phi[i] - 0.5 * sigmar * deltal
      end
    end
  end

  return
end

function apply_3d_muscl!(
  phi::OffsetArray{<:AbstractFloat, 3},
  phitilde::OffsetArray{<:AbstractFloat, 5},
  sizex::Integer,
  sizey::Integer,
  sizez::Integer,
  limitertype::MCVariant,
)
  phix = zeros(sizex)
  phiy = zeros(sizey)
  phiz = zeros(sizez)

  phitildex = zeros(sizex, 2)
  phitildey = zeros(sizey, 2)
  phitildez = zeros(sizez, 2)

  for kz in 2:(sizez - 1)
    for jy in 2:(sizey - 1)
      phix[:] .= phi.parent[:, jy, kz]
      muscle_reconstruct1D_mcvariant!(phix, phitildex, sizex)
      phitilde.parent[:, jy, kz, 1, :] .= phitildex.parent
    end
  end

  for kz in 2:(sizez - 1)
    for ix in 2:(sizex - 1)
      phiy .= phi.parent[ix, :, kz]
      muscle_reconstruct1D_mcvariant!(phiY, phitildey, sizey)
      phitilde.parent[ix, :, kz, 2, :] .= phitildey.parent
    end
  end

  for jy in 2:(sizey - 1)
    for ix in 2:(sizex - 1)
      phiz .= phi.parent[ix, jy, :]
      muscle_reconstruct1D_mcvariant!(phiz, phitildez, sizez)
      phitilde.parent[ix, jy, :, 3, :] .= phitildez.parent
    end
  end

  return
end

function reconstruct_rho!(state::State)
  (; sizex, sizey, sizez, nbx, nby, nbz) = state.namelists.domain
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
  apply_3d_muscl!(rhobar, rhotilde, sizex, sizey, sizez, limitertype)

  return
end

function reconstruct_rhop!(state::State)
  (; sizex, sizey, sizez, nbx, nby, nbz) = state.namelists.domain
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
  apply_3d_muscl!(rhopbar, rhoptilde, sizex, sizey, sizez, limitertype)

  return
end

function reconstruct_u!(state::State)
  (; sizex, sizey, sizez, nbx, nby, nbz) = state.namelists.domain
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

  apply_3d_muscl!(ubar, utilde, sizex, sizey, sizez, limitertype)

  return
end

function reconstruct_v!(state::State)
  (; sizex, sizey, sizez, nbx, nby, nbz) = state.namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = state.domain
  (; rho, v) = state.variables.predictands
  (; vbar) = state.variables.auxiliaries
  (; vtilde) = state.variables.reconstructions
  (; rhostrattfc, pstrattfc) = state.atmosphere

  for ix in (-nbx):(nx + nbx - 1)
    for jy in (-nby):(ny + nby)
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

  apply_3d_muscl!(vbar, vtilde, sizex, sizey, sizez, limitertype)

  return
end

function reconstruct_w!(state::State)
  (; namelists, domain) = state
  (; sizex, sizey, sizez, nbx, nby, nbz, nprocx, nprocy) = namelists.domain
  (; limitertype) = state.namelists.discretization
  (; nx, ny, nz) = domain
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
        wbar[ix, jy, kz] = wbar[ix, jy, kz] * rhoedgeu / pedgeu
      end
    end
  end

  apply_3d_muscl!(wbar, wtilde, sizex, sizey, sizez, limitertype)

  return
end

function reconstruct!(state::State)
  reconstruct_rho!(state)
  reconstruct_rhop!(state)
  reconstruct_u!(state)
  reconstruct_v!(state)
  reconstruct_w!(state)
  return
end
