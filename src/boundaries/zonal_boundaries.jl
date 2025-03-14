function set_zonal_boundaries_of_field!(
  field::OffsetMatrix{<:AbstractFloat},
  namelists::Namelists,
  domain::Domain,
)

  (; nbx, nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      field[nx + i, :] = field[i, :]
      field[-i + 1, :] = field[nx - j + 1, :]
    end
  end

  return
end

function set_zonal_boundaries_of_field!(
  field::OffsetArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)

  (; nbx, nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      field[nx + i, :, :] = field[i, :, :]
      field[-i + 1, :, :] = field[nx - j + 1, :, :]
    end
  end

  return
end

function set_zonal_boundaries_of_field!(
  field::OffsetArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)

  (; nbx, nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      field[nx + i, :, :, :] = field[i, :, :, :]
      field[-i + 1, :, :, :] = field[nx - j + 1, :, :, :]
    end
  end

  return
end

function set_zonal_boundaries!(
  state::State,
  variables::BoundaryPredictands,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (; rho, rhop, u, v, w, pip) = state.variables.predictands

  set_zonal_boundaries_of_field!(rho, namelists, domain)
  set_zonal_boundaries_of_field!(rhop, namelists, domain)
  set_zonal_boundaries_of_field!(u, namelists, domain)
  set_zonal_boundaries_of_field!(v, namelists, domain)
  set_zonal_boundaries_of_field!(w, namelists, domain)
  set_zonal_boundaries_of_field!(pip, namelists, domain)

  return
end

function set_zonal_boundaries!(
  state::State,
  variables::BoundaryReconstructions,
  boundaries::PeriodicBoundaries,
)
  (; namelists, domain) = state
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) =
    state.variables.reconstructions

  set_zonal_boundaries_of_field!(rhotilde, namelists, domain)
  set_zonal_boundaries_of_field!(rhoptilde, namelists, domain)
  set_zonal_boundaries_of_field!(utilde, namelists, domain)
  set_zonal_boundaries_of_field!(vtilde, namelists, domain)
  set_zonal_boundaries_of_field!(wtilde, namelists, domain)

  return
end
