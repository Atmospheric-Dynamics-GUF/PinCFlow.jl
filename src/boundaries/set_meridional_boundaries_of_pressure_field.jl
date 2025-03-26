function set_meridional_boundaries_of_pressure_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nprocy) = namelists.domain
  (; ny) = domain

  if nprocy > 1
    set_meridional_halos_of_pressure_field!(field, domain)
  else
    @views field[:, ny + 1, :] .= field[:, 1, :]
    @views field[:, 0, :] .= field[:, ny, :]
  end

  return
end
