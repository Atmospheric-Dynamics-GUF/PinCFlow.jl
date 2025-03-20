function set_meridional_boundaries_of_field!(
  field::AbstractMatrix{<:AbstractFloat},
  namelists::Namelists,
  domain::Domain,
)
  (; nby, nprocy) = namelists.domain
  (; ny) = domain
  if nprocy > 1
    set_meridional_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      field[:, ny + j] = view(field, :, j)
      field[:, -j + 1] = view(field, :, ny - j + 1)
    end
  end

  return
end

function set_meridional_boundaries_of_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nby, nprocy) = namelists.domain
  (; ny) = domain

  if nprocy > 1
    set_meridional_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      field[:, ny + j, :] .= view(field, :, j, :)
      field[:, -j + 1, :] .= view(field, :, ny - j + 1, :)
    end
  end

  return
end

function set_meridional_boundaries_of_field!(
  field::AbstractArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)
  (; nby, nprocy) = namelists.domain
  (; ny) = domain

  if nprocy > 1
    set_meridional_halos_of_field!(field, namelists, domain)
  else
    @inbounds for j in 1:nby
      field[:, ny + j, :, :, :] .= view(field, :, j, :, :, :)
      field[:, -j + 1, :, :, :] .= view(field, :, ny - j + 1, :, :, :)
    end
  end

  return
end
