function set_zonal_boundaries!(
  state::State,
  variables::BoundaryPredictands,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (; nbx, nprocx) = namelists.domain
  (; nx) = domain
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  if nprocx > 1
    set_zonal_halos_of_field!(rho, namelists, domain)
    set_zonal_halos_of_field!(rhop, namelists, domain)
    set_zonal_halos_of_field!(u, namelists, domain)
    set_zonal_halos_of_field!(v, namelists, domain)
    set_zonal_halos_of_field!(w, namelists, domain)
    set_zonal_halos_of_field!(pip, namelists, domain)
  else
    for i in 1:nbx
      rho[nx + i, :, :] = rho[i, :, :]
      rho[-i + 1, :, :] = rho[nx - i + 1, :, :]

      rhop[nx + i, :, :] = rhop[i, :, :]
      rhop[-i + 1, :, :] = rhop[nx - i + 1, :, :]

      u[nx + i, :, :] = u[i, :, :]
      u[-i + 1, :, :] = u[nx - i + 1, :, :]

      v[nx + i, :, :] = v[i, :, :]
      v[-i + 1, :, :] = v[nx - i + 1, :, :]

      w[nx + i, :, :] = w[i, :, :]
      w[-i + 1, :, :] = w[nx - i + 1, :, :]

      pip[nx + i, :, :] = pip[i, :, :]
      pip[-i + 1, :, :] = pip[nx - i + 1, :, :]
    end
  end

  return
end

function set_zonal_boundaries!(
  state::State,
  variables::BoundaryReconstructions,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (; nbx, nprocx) = namelists.domain
  (; nx) = domain
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) =
    state.variables.reconstructions

  if nprocx > 1
    set_zonal_halos_of_field!(rhotilde, namelists, domain)
    set_zonal_halos_of_field!(rhoptilde, namelists, domain)
    set_zonal_halos_of_field!(utilde, namelists, domain)
    set_zonal_halos_of_field!(vtilde, namelists, domain)
    set_zonal_halos_of_field!(wtilde, namelists, domain)
  else
    for i in 1:nbx
      rhotilde[nx + i, :, :, :, :] = rhotilde[i, :, :, :, :]
      rhotilde[-i + 1, :, :] = rhotilde[nx - i + 1, :, :, :, :]

      rhoptilde[nx + i, :, :, :, :] = rhoptilde[i, :, :, :, :]
      rhoptilde[-i + 1, :, :, :, :] = rhoptilde[nx - i + 1, :, :, :, :]

      utilde[nx + i, :, :, :, :] = utilde[i, :, :, :, :]
      utilde[-i + 1, :, :, :, :] = utilde[nx - i + 1, :, :, :, :]

      vtilde[nx + i, :, :, :, :] = vtilde[i, :, :, :, :]
      vtilde[-i + 1, :, :, :, :] = vtilde[nx - i + 1, :, :, :, :]

      wtilde[nx + i, :, :, :, :] = wtilde[i, :, :, :, :]
      wtilde[-i + 1, :, :, :, :] = wtilde[nx - i + 1, :, :, :, :]
    end
  end

  return
end
