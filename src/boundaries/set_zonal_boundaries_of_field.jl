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
