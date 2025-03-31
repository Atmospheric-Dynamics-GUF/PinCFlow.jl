function set_meridional_boundaries_of_pressure_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nprocy) = namelists.domain
  (; i0, i1, j0, j1, k0, k1) = domain

  if nprocy > 1
    set_meridional_halos_of_pressure_field!(field, domain)
  else
    @views field[(i0 - 1):(i1 + 1), j0 - 1, (k0 - 1):(k1 + 1)] .=
      field[(i0 - 1):(i1 + 1), j1, (k0 - 1):(k1 + 1)]
    @views field[(i0 - 1):(i1 + 1), j1 + 1, (k0 - 1):(k1 + 1)] .=
      field[(i0 - 1):(i1 + 1), j0, (k0 - 1):(k1 + 1)]
  end

  return
end
