function remove_ray_volumes!(state::State)
  (; i0, i1, j0, j1) = state.domain
  (; nray, rays) = state.wkb
  for kz in k0:k1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
    if nray[ix, jy, kz] <= 0
      continue
    end
    nrlc = 0
    for iray in 1:nray[ix, jy, kz]
      if rays.dens[iray, ix, jy, kz] == 0.0
        continue
      end
      nrlc = nrlc + 1
      copy_ray!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))
    end
    nray[ix, jy, kz] = nrlc
  end
  return
end