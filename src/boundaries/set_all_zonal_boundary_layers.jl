function set_all_zonal_boundary_layers!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)
  (; nbx, nprocx) = namelists.domain
  (; i0, i1) = domain

  if nprocx > 1
    set_all_zonal_halo_layers!(field, namelists, domain)
  else
    for i in 1:nbx
      @views field[i0 - i, :, :] .= field[i1 - i + 1, :, :]
      @views field[i1 + i, :, :] .= field[i0 + i - 1, :, :]
    end
  end

  return
end

function set_all_zonal_boundary_layers!(
  field::AbstractArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)
  (; nbx, nprocx) = namelists.domain
  (; i0, i1) = domain

  if nprocx > 1
    set_all_zonal_halo_layers!(field, namelists, domain)
  else
    for i in 1:nbx
      @views field[i0 - i, :, :, :, :] .= field[i1 - i + 1, :, :, :, :]
      @views field[i1 + i, :, :, :, :] .= field[i0 + i - 1, :, :, :, :]
    end
  end

  return
end
