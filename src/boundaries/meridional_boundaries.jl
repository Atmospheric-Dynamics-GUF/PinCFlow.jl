function set_meridional_boundaries!(
  state::State,
  variables::BoundaryPredictands,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (; nby, nprocy) = namelists.domain
  (; ny) = domain
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  if nprocy > 1
    set_meridional_halos_of_field!(rho, namelists, domain)
    set_meridional_halos_of_field!(rhop, namelists, domain)
    set_meridional_halos_of_field!(u, namelists, domain)
    set_meridional_halos_of_field!(v, namelists, domain)
    set_meridional_halos_of_field!(w, namelists, domain)
    set_meridional_halos_of_field!(pip, namelists, domain)
  else
    for j in 1:nby
      rho[:, ny + j, :] = rho[:, j, :]
      rho[:, -j + 1, :] = rho[:, ny - j + 1, :]

      rhop[:, ny + j, :] = rhop[:, j, :]
      rhop[:, -j + 1, :] = rhop[:, ny - j + 1, :]

      u[:, ny + j, :] = u[:, j, :]
      u[:, -j + 1, :] = u[:, ny - j + 1, :]

      v[:, ny + j, :] = v[:, j, :]
      v[:, -j + 1, :] = v[:, ny - j + 1, :]

      w[:, ny + j, :] = w[:, j, :]
      w[:, -j + 1, :] = w[:, ny - j + 1, :]

      pip[:, ny + j, :] = pip[:, j, :]
      pip[:, -j + 1, :] = pip[:, ny - j + 1, :]
    end
  end

  return
end

function set_meridional_boundaries!(
  state::State,
  variables::BoundaryReconstructions,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (; nby, nprocy) = namelists.domain
  (; ny) = domain
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) =
    state.variables.reconstructions

  if nprocy > 1
    set_meridional_halos_of_field!(rhotilde, namelists, domain)
    set_meridional_halos_of_field!(rhoptilde, namelists, domain)
    set_meridional_halos_of_field!(utilde, namelists, domain)
    set_meridional_halos_of_field!(vtilde, namelists, domain)
    set_meridional_halos_of_field!(wtilde, namelists, domain)
  else
    for i in 1:nbx
      rhotilde[:, ny + j, :, :, :] = rhotilde[:, j, :, :, :]
      rhotilde[:, -j + 1, :, :, :] = rhotilde[:, ny - j + 1, :, :, :]

      rhoptilde[:, ny + j, :, :, :] = rhoptilde[:, j, :, :, :]
      rhoptilde[:, -j + 1, :, :, :] = rhoptilde[:, ny - j + 1, :, :, :]

      utilde[:, ny + j, :, :, :] = utilde[:, j, :, :, :]
      utilde[:, -j + 1, :, :, :] = utilde[:, ny - j + 1, :, :, :]

      vtilde[:, ny + j, :, :, :] = vtilde[:, j, :, :, :]
      vtilde[:, -j + 1, :, :, :] = vtilde[:, ny - j + 1, :, :, :]

      wtilde[:, ny + j, :, :, :] = wtilde[:, j, :, :, :]
      wtilde[:, -j + 1, :, :, :] = wtilde[:, ny - j + 1, :, :, :]
    end
  end

  return
end
