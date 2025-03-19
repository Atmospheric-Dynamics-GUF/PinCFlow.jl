function set_zonal_boundaries_of_field!(
  field::AbstractMatrix{<:AbstractFloat},
  namelists::Namelists,
  domain::Domain,
)
  (; nbx, nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_field!(field, namelists, domain)
  else
    for i in 1:nbx
      field[nx + i, :] = field[i, :]
      field[-i + 1, :] = field[nx - i + 1, :]
    end
  end

  return
end

function set_zonal_boundaries_of_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nbx, nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_field!(field, namelists, domain)
  else
    for i in 1:nbx
      field[nx + i, :, :] = field[i, :, :]
      field[-i + 1, :, :] = field[nx - i + 1, :, :]
    end
  end

  return
end

function set_zonal_boundaries_of_field!(
  field::AbstractArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)
  (; nbx, nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_field!(field, namelists, domain)
  else
    for i in 1:nbx
      field[nx + i, :, :, :, :] = field[i, :, :, :, :]
      field[-i + 1, :, :, :, :] = field[nx - i + 1, :, :, :, :]
    end
  end

  return
end
