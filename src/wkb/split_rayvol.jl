function split_z(ix, jy, kz, rays, nray, dz)
  nrlc = nray[ix, jy, kz]
  for iray in 1:nray[ix, jy, kz]
      ijk = (iray, ix, jy, kz)

      zr = rays.z[(ijk)...]

    dzr = ray.dzray[(ijk)...]
    azm = ray.area_zm[(ijk)...]

    xr = ray.x[(ijk)...]
    yr = ray.y[(ijk)...]

    ixrv = nint((xr - lx[0]) / dx + 0.5) - ix0
    jyrv = nint((yr - ly[0]) / dy + 0.5) - jy0
    kzrvd = kzTildeTFC[ixrv, jyrv, zr - 0.5 * dzr]
    kzrvu = kzTildeTFC[ixrv, jyrv, zr + 0.5 * dzr]

    dzmin = dz
    for kzrv in kzrvd:kzrvu
      dzmin = min(dzmin, jac[ixrv, jyrv, kzrv] * dz)
    end

    if (dzr > dzmin)
      factor = ceil(dzr / dzmin)
      rays.z[(ijk)...] = zr + 0.5 * (1 / factor - 1) * dzr
      rays.dzray[(ijk)...] = dzr / factor
      rays.area_zm[(ijk)...] = azm / factor
      for jRay in (nrlc + 1):(nrlc + factor - 1)
        ray[jRay, ix, jy, kz] = ray[iRay, ix, jy, kz]
        rays.z[(ijk)...] += (jRray - nrlc) * dzr / factor
      end
      nrlc = nrlc + factor - 1
    end
  end

  if (nrlc > nray[ix, jy, kz])
    nray[ix, jy, kz] = nrlc

    if (nray[ix, jy, kz] > nray_wrk)
      # print *, 'ERROR at ix,jy,kz =', ix, jy, kz
      # print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
      exit()
    end
  end
end

function split_y(ix, jy, kz, rays, nray, dy)
  nrlc = nray[ix, jy, kz]
  for iray in 1:nray[ix, jy, kz]
    ray = rays(iray, ix, jy, kz)
    yr = ray.y
    dyr = ray.dyray

    if (dyr > dy)
      nrlc = nrlc + 1

      if (nrlc > nray_wrk)
        # print *, 'r.v. getting too many by splitting in y  direction'
        exit()
      end

      ayl = ray.area_yl

      ray.dyray = 0.5 * dyr
      ray.area_yl = 0.5 * ayl

      # copy ?
      rays[nrlc, ix, jy, kz] = ray

      ray.y = yr - 0.25 * dyr
      ray[nrlc, ix, jy, kz].y = yr + 0.25 * dyr
    end
  end

  if (nrlc > nray[ix, jy, kz])
    nray[ix, jy, kz] = nrlc
  end
    if (nray[ix, jy, kz] > nray_wrk)
      # print *, 'ERROR at ix,jy,kz =', ix, jy, kz
      # print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
      exit()
    end
  end

function split_x(ix, jy, kz, rays, nray, dx)
  nrlc = nray[ix, jy, kz]
  for iRay in 1:nray[ix, jy, kz]
    ray = rays[iray, ix, jy, kz]
    dxr = ray.dxray

    if (dxr > dx)
      nrlc = nrlc + 1

      if (nrlc > nray_wrk)
        # print *, 'r.v. getting too many by splitting in x  direction'
        exit()
      end

      xr = ray.x

      axk = ray.area_xk

      ray.dxray = 0.5 * dxr

      ray.area_xk = 0.5 * axk

      rays[nrlc, ix, jy, kz] = ray

      ray.x = xr - 0.25 * dxr
      rays[nrlc, ix, jy, kz].x = xr + 0.25 * dxr
    end
  end
  if (nrlc > nray[ix, jy, kz])
    nray[ix, jy, kz] = nrlc

    if (nRay[ix, jy, kz] > nray_wrk)
      # print *, 'ERROR at ix,jy,kz =', ix, jy, kz
      # print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
      exit()
    end
  end
end

function split_rayvol(state)
  if (state.steady_state)
    return
  end

  ix0 = is + nbx - 1
  jy0 = js + nby - 1

  # ! total number of ray volumes before splitting

  # nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

  # ! testb
  # ! print*,'before splitting nrvloc =',nrvloc
  # ! teste

  # call mpi_reduce(nrvloc, nrvtt0, 1, mpi_integer, mpi_sum, root, comm, ierror)
  # call mpi_bcast(nrvtt0, 1, mpi_integer, root, comm, ierror)
  (; nx, ny, nz) = state.domain
  domain = state.domain
  (; dx, dy, dz) = state.grid
  (; rays, nray) = state.wkb
  for kz in 1:nz
    for jy in 1:ny
      for ix in 1:nx
        if (nray[ix, jy, kz] < 1)
          continue
        end

        if (domain.sizex > 1)
          split_x(ix, jy, kz, rays, nray, dx)
        end
        if (domain.sizey > 1)
          split_y(ix, jy, kz, rays, nray, dy)
        end
        split_z(ix, jy, kz, rays, nray, dz)
      end
    end
  end

  # ! total number of ray volumes after splitting

  # nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

  # ! testb
  # ! print*,'after splitting nrvloc =',nrvloc
  # ! teste

  # call mpi_reduce(nrvloc, nrvtt1, 1, mpi_integer, mpi_sum, root, comm, ierror)
  # call mpi_bcast(nrvtt1, 1, mpi_integer, root, comm, ierror)

  # if(master .and. nrvtt1 > nrvtt0) then
  #   print *, 'after splitting nray =', nrvtt1
  # end if

  return
end
