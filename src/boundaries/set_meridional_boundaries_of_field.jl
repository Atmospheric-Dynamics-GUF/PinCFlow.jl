function set_meridional_boundaries_of_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nby, nprocy) = namelists.domain
  (; j0, j1) = domain

  if nprocy > 1
    set_meridional_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      @views field[:, j0 - j, :] .= field[:, j1 - j + 1, :]
      @views field[:, j1 + j, :] .= field[:, j0 + j - 1, :]
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
  (; j0, j1) = domain

  if nprocy > 1
    set_meridional_halos_of_field!(field, namelists, domain)
  else
    for j in 1:nby
      @views field[:, j0 - j, :, :, :] .= field[:, j1 - j + 1, :, :, :]
      @views field[:, j1 + j, :, :, :] .= field[:, j0 + j - 1, :, :, :]
    end
  end

  return
end
