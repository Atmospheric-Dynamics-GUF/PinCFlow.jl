function set_one_zonal_boundary_layer!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nprocx) = namelists.domain
  (; i0, i1, j0, j1, k0, k1) = domain

  if nprocx > 1
    set_one_zonal_halo_layer!(field, domain)
  else
    @views field[i0 - 1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .=
      field[i1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]
    @views field[i1 + 1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .=
      field[i0, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]
  end

  return
end
