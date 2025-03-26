function set_zonal_boundaries_of_pressure_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nprocx) = namelists.domain
  (; nx) = domain

  if nprocx > 1
    set_zonal_halos_of_pressure_field!(field, domain)
  else
    @views field[nx + 1, :, :] .= field[1, :, :]
    @views field[0, :, :] .= field[nx, :, :]
  end

  return
end
