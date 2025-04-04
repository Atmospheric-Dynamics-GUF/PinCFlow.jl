struct X end
struct Y end
struct Z end

function shift_ray_volumes!(state::State)
  (; sizex, sizey, nprocx, nprocy) = state.namelists.domain
  (; zboundaries) = state.namelists.setting
  (; steady_state) = state.namelists.wkb

  if steady_state
    return
  else
    if sizex > 1
      if nprocx > 1
        error("Parallelized shfiting not ready yet!")
      else
        set_zonal_boundary_ray_volumes!(state)
        shift_ray_volumes!(state, X())
        set_zonal_boundary_ray_volumes!(state)
        remove_ray_volumes!(state)
      end
    end

    if sizey > 1
      if nprocy > 1
        error("Parallelized shfiting not ready yet!")
      else
        set_meridional_boundary_ray_volumes!(state)
        shift_ray_volumes!(state, Y())
        set_meridional_boundary_ray_volumes!(state)
        remove_ray_volumes!(state)
      end
    end

    set_vertical_boundary_ray_volumes(state, zboundaries)
    shift_ray_volumes!(state, Z())
    set_vertical_boundary_ray_volumes(state, zboundaries)
    remove_ray_volumes!(state)

    check_ray_volumes(state)
    return
  end
end

function shift_ray_volumes!(state::State, direction::X)
  (; io, i0, i1, j0, j1, k0, k1) = state.domain
  (; lx, dx) = state.grid
  (; nray, rays) = state.wkb

  for kzrv in (k0 - 1):(k1 + 1), jyrv in (j0 - 1):(j1 + 1), ixrv in i0:i1
    for iray in 1:nray[ixrv, jyrv, kzrv]
      if rays.dens[iray, ixrv, jyrv, kzrv] != 0.0
        xr = rays.x[iray, ixrv, jyrv, kzrv]
        ix = round(Int, (xr - lx[1] - dx / 2) / dx) + i0 - io

        if ix != ixrv
          nray[ixrv, jyrv, kz] += 1
          jray = nray[ixrv, jyrv, kz]
          if jray > nray_wrk
            error("Error in shift_ray_volumes!: nray > nray_wrk!")
          end
          copy_ray!(rays, (iray, ixrv, jyrv, kzrv), (jray, ixrv, jyrv, kz))
          rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
        end
      end
    end
  end
end

function shift_ray_volumes!(state::State, direction::Y)
  (; jo, i0, i1, j0, j1, k0, k1) = state.domain
  (; ly, dy) = state.grid
  (; nray, rays) = state.wkb

  for kzrv in (k0 - 1):(k1 + 1), jyrv in j0:j1, ixrv in (i0 - 1):(i1 + 1)
    for iray in 1:nray[ixrv, jyrv, kzrv]
      if rays.dens[iray, ixrv, jyrv, kzrv] != 0.0
        yr = rays.y[iray, ixrv, jyrv, kzrv]
        jy = round(Int, (yr - ly[1] - dy / 2) / dy) + j0 - jo

        if jy != jyrv
          nray[ixrv, jyrv, kz] += 1
          jray = nray[ixrv, jyrv, kz]
          if jray > nray_wrk
            error("Error in shift_ray_volumes!: nray > nray_wrk!")
          end
          copy_ray!(rays, (iray, ixrv, jyrv, kzrv), (jray, ixrv, jyrv, kz))
          rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
        end
      end
    end
  end
end

function shift_ray_volumes!(state::State, direction::Z)
  (; domain, grid) = state
  (; i0, i1, j0, j1, k0, k1) = domain
  (; nray, rays) = state.wkb

  for kzrv in k0:k1, jyrv in (j0 - 1):(j1 + 1), ixrv in (i0 - 1):(i1 + 1)
    for iray in 1:nray[ixrv, jyrv, kzrv]
      if rays.dens[iray, ixrv, jyrv, kzrv] != 0.0
        zr = rays.z[iray, ixrv, jyrv, kzrv]
        kz = kztildetfc(ixrv, jyrv, zr, domain, grid)

        if kz > k1
          kz = k1
        end
        if kz < k0
          kz = k0
        end

        if kz != kzrv
          nray[ixrv, jyrv, kz] += 1
          jray = nray[ixrv, jyrv, kz]
          if jray > nray_wrk
            error("Error in shift_ray_volumes!: nray > nray_wrk!")
          end
          copy_ray!(rays, (iray, ixrv, jyrv, kzrv), (jray, ixrv, jyrv, kz))
          rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
        end
      end
    end
  end

  return
end
