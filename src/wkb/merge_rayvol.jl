
function wave_numbers_interval(ix, jy, kz, rays, nray)
  for iRay in 1:nRay[ix, jy, kz]
    ray = rays[iRai, ik, jy, kz]

    wnrk = ray.k
    dwnrk = abs(ray.dkray)

    wnrl = ray.l
    dwnrl = abs(ray.dlray)

    wnrm = ray.m
    dwnrm = abs(ray.dmray)

    if (sizeX > 1)

      # ! subdivide the interval of wavenumbers to be processed
      # ! into up to three ranges:
      # ! (1) if there are negative wavenumbers wnrk < 0, an
      # ! interval -wnrk_max_n <= wnrk <= -wnrk_min_n
      # ! (2) if there are zero wavnumbers wnrk = 0, a slot
      # ! for those
      # ! (3) if there are positive wavenumbers wnrk > 0, an
      # ! interval wnrk_min_p <= wnrk <= wnrk_max_p

      if (wnrk > 0.0)
        if (wnrk_min_p == 0.0)
          wnrk_min_p = wnrk
        else
          wnrk_min_p = min(wnrk_min_p, wnrk)
        end

        wnrk_max_p = max(wnrk_max_p, wnrk)
      elseif (wnrk < 0.0)
        if (wnrk_min_n == 0.0)
          wnrk_min_n = -wnrk
        else
          wnrk_min_n = min(wnrk_min_n, -wnrk)
        end

        wnrk_max_n = max(wnrk_max_n, -wnrk)
      end
    end

    if (sizeY > 1)
      if (wnrl > 0.0)
        if (wnrl_min_p == 0.0)
          wnrl_min_p = wnrl
        else
          wnrl_min_p = min(wnrl_min_p, wnrl)
        end

        wnrl_max_p = max(wnrl_max_p, wnrl)
      elseif (wnrl < 0.0)
        if (wnrl_min_n == 0.0)
          wnrl_min_n = -wnrl
        else
          wnrl_min_n = min(wnrl_min_n, -wnrl)
        end

        wnrl_max_n = max(wnrl_max_n, -wnrl)
      end
    end

    if (wnrm > 0.0)
      if (wnrm_min_p == 0.0)
        wnrm_min_p = wnrm
      else
        wnrm_min_p = min(wnrm_min_p, wnrm)
      end

      wnrm_max_p = max(wnrm_max_p, wnrm)
    elseif (wnrm < 0.0)
      if (wnrm_min_n == 0.0)
        wnrm_min_n = -wnrm
      else
        wnrm_min_n = min(wnrm_min_n, -wnrm)
      end

      wnrm_max_n = max(wnrm_max_n, -wnrm)
    end
  end

  if (sizeX > 1)
    if (wnrk_min_n == 0.0 && wnrk_max_n == 0.0)
      if (wnrk_min_p /= 0.0 && wnrk_max_p /= 0.0)
        wnrk_min_n = wnrk_min_p
        wnrk_max_n = wnrk_max_p
      else
        # ! all limits zero only applies if all wnrk = 0
        # ! hence, just in order to provide some numbers ...
        wnrk_min_n = 1.0
        wnrk_max_n = 2.0
      end
    end

    if (wnrk_min_p == 0.0 && wnrk_max_p == 0.0)
      if (wnrk_min_n /= 0.0 && wnrk_max_n /= 0.0)
        wnrk_min_p = wnrk_min_n
        wnrk_max_p = wnrk_max_n
      else
        # ! all limits zero only applies if all wnrk = 0
        # ! hence, just in order to provide some numbers ...
        wnrk_min_p = 1.0
        wnrk_max_p = 2.0
      end
    end

    # ! in order to prevent zero-width intervals ...

    if (wnrk_min_n == wnrk_max_n)
      wnrk_min_n = 0.5 * wnrk_min_n
      wnrk_max_n = 2.0 * wnrk_max_n
    end

    if (wnrk_min_p == wnrk_max_p)
      wnrk_min_p = 0.5 * wnrk_min_p
      wnrk_max_p = 2.0 * wnrk_max_p
    end
  end
  if (sizeY > 1)
    if (wnrl_min_n == 0.0 && wnrl_max_n == 0.0)
      if (wnrl_min_p /= 0.0 && wnrl_max_p /= 0.0)
        wnrl_min_n = wnrl_min_p
        wnrl_max_n = wnrl_max_p
      else
        # ! all limits zero only applies if all wnrl = 0
        # ! hence, just in order to provide some numbers ...
        wnrl_min_n = 1.0
        wnrl_max_n = 2.0
      end
    end

    if (wnrl_min_p == 0.0 && wnrl_max_p == 0.0)
      if (wnrl_min_n /= 0.0 && wnrl_max_n /= 0.0)
        wnrl_min_p = wnrl_min_n
        wnrl_max_p = wnrl_max_n
      else
        # ! all limits zero only applies if all wnrl = 0
        # ! hence, just in order to provide some numbers ...
        wnrl_min_p = 1.0
        wnrl_max_p = 2.0
      end
    end

    # ! in order to prevent zero-width intervals ...

    if (wnrl_min_n == wnrl_max_n)
      wnrl_min_n = 0.5 * wnrl_min_n
      wnrl_max_n = 2.0 * wnrl_max_n
    end

    if (wnrl_min_p == wnrl_max_p)
      wnrl_min_p = 0.5 * wnrl_min_p
      wnrl_max_p = 2.0 * wnrl_max_p
    end
  end # sizey > 1

  # z
  if (wnrm_min_n == 0.0 && wnrm_max_n == 0.0)
    if (wnrm_min_p /= 0.0 && wnrm_max_p /= 0.0)
      wnrm_min_n = wnrm_min_p
      wnrm_max_n = wnrm_max_p
    else
      # ! all limits zero only applies if all wnrm = 0
      # ! hence, just in order to provide some numbers ...
      wnrm_min_n = 1.0
      wnrm_max_n = 2.0
    end
  end

  if (wnrm_min_p == 0.0 && wnrm_max_p == 0.0)
    if (wnrm_min_n /= 0.0 && wnrm_max_n /= 0.0)
      wnrm_min_p = wnrm_min_n
      wnrm_max_p = wnrm_max_n
    else
      # ! all limits zero only applies if all wnrm = 0
      # ! hence, just in order to provide some numbers ...
      wnrm_min_p = 1.0
      wnrm_max_p = 2.0
    end
  end

  # ! in order to prevent zero-width intervals ...

  if (wnrm_min_n == wnrm_max_n)
    wnrm_min_n = 0.5 * wnrm_min_n
    wnrm_max_n = 2.0 * wnrm_max_n
  end

  if (wnrm_min_p == wnrm_max_p)
    wnrm_min_p = 0.5 * wnrm_min_p
    wnrm_max_p = 2.0 * wnrm_max_p
  end
  return (;)
end

function log_widths(max_n, min_n, max_p, min_p, nray)
  mg_n = log(max_n / min_n) / (nray / 2 - 1)
  mg_p = log(max_p / min_p) / (nray / 2 - 1)
  return (mg_n, mg_p)
end

function generate_merged_rv(rays, nRray)
  for iRay in 1:nRay(ix, jy, kz)
    ray = rays[iRay, ix, jy, kz]
    wnrk = ray.k
    dwnrk = ray.dkray

    wnrl = ray.l
    dwnrl = ray.dlray

    wnrm = ray.m
    dwnrm = ray.dmray

    xr = ray.x
    dxr = ray.dxray

    yr = ray.y
    dyr = ray.dyray

    zr = ray.z
    dzr = ray.dzray

    axk = ray.area_xk
    ayl = ray.area_yl
    azm = ray.area_zm

    wdr = ray.dens
    omir = ray.omega
    if (sizeX > 1)
      fcpspx = axk

      # ! nxRay - 1 intervals for k:
      # ! indices 1 ... nxRay/2 - 1 for negative k
      # ! index nxRay/2 for k = 0
      # ! indices nxRay/2 + 1 ... nxRay - 1 for positive k

      if (wnrk < 0.0)
        if (abs(log(-wnrk / wnrk_max_n) / dwnrk_mg_n) < 1.0)
          ir_k = nxRay / 2 - 1
        else
          ir_k = int(log(-wnrk / wnrk_min_n) / dwnrk_mg_n) + 1
        end
      elseif (wnrk == 0.0)
        ir_k = nxRay / 2
      else
        if (abs(log(wnrk / wnrk_max_p) / dwnrk_mg_p) < 1.0)
          ir_k = nxRay - 1
        else
          ir_k = int(log(wnrk / wnrk_min_p) / dwnrk_mg_p) + nxRay / 2 + 1
        end
      end

      if (ir_k < 1)
        # print *, 'ERROR in merge_rayvol: ir_k =', ir_k, '< 1'
        # print *, 'wnrk = ', wnrk
        # print *, 'wnrk_min_n = ', wnrk_min_n
        # print *, 'wnrk_max_n = ', wnrk_max_n
        # print *, 'wnrk_min_p = ', wnrk_min_p
        # print *, 'wnrk_max_p = ', wnrk_max_p
        # stop
        exit()
      elseif (ir_k > nxRay - 1)
        # print *, 'ERROR in merge_rayvol: ir_k =', ir_k, '> nxRay - 1 &
        #     &=', nxRay - 1
        # print *, 'wnrk = ', wnrk
        # print *, 'wnrk_min_n = ', wnrk_min_n
        # print *, 'wnrk_max_n = ', wnrk_max_n
        # print *, 'wnrk_min_p = ', wnrk_min_p
        # print *, 'wnrk_max_p = ', wnrk_max_p
        # stop
        exit()
      end
    else
      fcpspx = 1.0
      ir_k = 1
    end # sizeX > 1

    if (sizeY > 1)
      fcpspy = ayl

      if (wnrl < 0.0)
        if (abs(log(-wnrl / wnrl_max_n) / dwnrl_mg_n) < 1.0)
          ir_l = nyRay / 2 - 1
        else
          ir_l = int(log(-wnrl / wnrl_min_n) / dwnrl_mg_n) + 1
        end
      elseif (wnrl == 0.0)
        ir_l = nyRay / 2
      else
        if (abs(log(wnrl / wnrl_max_p) / dwnrl_mg_p) < 1.0)
          ir_l = nyRay - 1
        else
          ir_l = int(log(wnrl / wnrl_min_p) / dwnrl_mg_p) + nyRay / 2 + 1
        end
      end

      if (ir_l < 1)
        # print *, 'ERROR in merge_rayvol: ir_l =', ir_l, '< 1'
        # print *, 'wnrl = ', wnrl
        # print *, 'wnrl_min_n = ', wnrl_min_n
        # print *, 'wnrl_max_n = ', wnrl_max_n
        # print *, 'wnrl_min_p = ', wnrl_min_p
        # print *, 'wnrl_max_p = ', wnrl_max_p
        exit()
      elseif (ir_l > nyRay - 1)
        # print *, 'ERROR in merge_rayvol: ir_l =', ir_l, '> nyRay - 1 &
        #     &=', nyRay - 1
        # print *, 'wnrl = ', wnrl
        # print *, 'wnrl_min_n = ', wnrl_min_n
        # print *, 'wnrl_max_n = ', wnrl_max_n
        # print *, 'wnrl_min_p = ', wnrl_min_p
        # print *, 'wnrl_max_p = ', wnrl_max_p
        # stop
        exit()
      end
    else
      fcpspy = 1.0
      ir_l = 1
    end # sizey > 1

    fcpspz = azm

    if (wnrm < 0.0)
      if (abs(log(-wnrm / wnrm_max_n) / dwnrm_mg_n) < 1.0)
        ir_m = nzRay / 2 - 1
      else
        ir_m = int(log(-wnrm / wnrm_min_n) / dwnrm_mg_n) + 1
      end
    elseif (wnrm == 0.0)
      ir_m = nzRay / 2
    else
      if (abs(log(wnrm / wnrm_max_p) / dwnrm_mg_p) < 1.0)
        ir_m = nzRay - 1
      else
        ir_m = int(log(wnrm / wnrm_min_p) / dwnrm_mg_p) + nzRay / 2 + 1
      end
    end

    if (ir_m < 1)
      # print *, 'ERROR in merge_rayvol: ir_m =', ir_m, '< 1'
      # print *, 'wnrm = ', wnrm
      # print *, 'wnrm_min_n = ', wnrm_min_n
      # print *, 'wnrm_max_n = ', wnrm_max_n
      # print *, 'wnrm_min_p = ', wnrm_min_p
      # print *, 'wnrm_max_p = ', wnrm_max_p
      # print *, 'wnrm_max_p - wnrm = ', wnrm_max_p - wnrm
      exit()
    elseif (ir_m > nzRay - 1)
      then
      # print *, 'ERROR in merge_rayvol: ir_m =', ir_m, '> nzRay - 1 =', &
      #     &nzRay - 1
      # print *, 'wnrm = ', wnrm
      # print *, 'wnrm_min_n = ', wnrm_min_n
      # print *, 'wnrm_max_n = ', wnrm_max_n
      # print *, 'wnrm_min_p = ', wnrm_min_p
      # print *, 'wnrm_max_p = ', wnrm_max_p
      # print *, 'wnrm_max_p - wnrm = ', wnrm_max_p - wnrm
      exit()
    end

    if (sizeX > 1)
      if (sizeY > 1)
        jRay =
          (ir_m - 1) * (nyRay - 1) * (nxRay - 1) +
          (ir_l - 1) * (nxRay - 1) +
          ir_k
      else
        jRay = (ir_m - 1) * (nxRay - 1) + ir_k
      end
    else
      if (sizeY > 1)
        jRay = (ir_m - 1) * (nyRay - 1) + ir_l
      else
        jRay = ir_m
      end
    end

    nr_merge[jRay] += +1

    if (nr_merge[jRay] == 1)
      xrmnmg[jRay] = xr - 0.5 * dxr
      xrmxmg[jRay] = xr + 0.5 * dxr

      yrmnmg[jRay] = yr - 0.5 * dyr
      yrmxmg[jRay] = yr + 0.5 * dyr

      zrmnmg[jRay] = zr - 0.5 * dzr
      zrmxmg[jRay] = zr + 0.5 * dzr

      krmnmg[jRay] = wnrk - 0.5 * dwnrk
      krmxmg[jRay] = wnrk + 0.5 * dwnrk

      lrmnmg[jRay] = wnrl - 0.5 * dwnrl
      lrmxmg[jRay] = wnrl + 0.5 * dwnrl

      mrmnmg[jRay] = wnrm - 0.5 * dwnrm
      mrmxmg[jRay] = wnrm + 0.5 * dwnrm

      if (cons_merge == "wa")
        # ! wave-action density after merging to be determined
        # ! such that the wave action remains the same,
        # ! hence ...

        function wadrmg[jray]
          return wdr * fcpspx * fcpspy * fcpspz
        end
      elseif (cons_merge == "en")
        # ! wave-action density after merging to be determined
        # ! such that the wave energy remains the same,
        # ! hence ...

        function wadrmg[jray]
          return wdr * omir * fcpspx * fcpspy * fcpspz
        end
      else
        # stop 'wrong cons_merge in merge_rayvol'
        exit()
      end
    else
      xrmnmg[jRay] = min(xrmnmg[jRay], xr - 0.5 * dxr)
      xrmxmg[jRay] = max(xrmxmg[jRay], xr + 0.5 * dxr)

      yrmnmg[jray] = min(yrmnmg[jray], yr - 0.5 * dyr)
      yrmxmg[jray] = max(yrmxmg[jray], yr + 0.5 * dyr)

      zrmnmg[jray] = min(zrmnmg[jray], zr - 0.5 * dzr)
      zrmxmg[jray] = max(zrmxmg[jray], zr + 0.5 * dzr)

      krmnmg[jray] = min(krmnmg[jray], wnrk - 0.5 * dwnrk)
      krmxmg[jray] = max(krmxmg[jray], wnrk + 0.5 * dwnrk)

      lrmnmg[jray] = min(lrmnmg[jray], wnrl - 0.5 * dwnrl)
      lrmxmg[jray] = max(lrmxmg[jray], wnrl + 0.5 * dwnrl)

      mrmnmg[jray] = min(mrmnmg[jray], wnrm - 0.5 * dwnrm)
      mrmxmg[jray] = max(mrmxmg[jray], wnrm + 0.5 * dwnrm)

      if (cons_merge == "wa")
        # ! wave-action density after merging to be determined
        # ! such that the wave action remains the same,
        # ! hence ...

        wadrmg[jRay] = wadrmg[jRay] + wdr * fcpspx * fcpspy * fcpspz
      elseif (cons_merge == "en")
        # ! wave-action density after merging to be determined
        # ! such that the wave energy remains the same,
        # ! hence ...

        wadrmg[jRay] += wdr * omir * fcpspx * fcpspy * fcpspz
      else
        # stop 'wrong cons_merge in merge_rayvol'
        exit()
      end
    end
  end
  return
end

function replace_rayvolume(nr_merge, rays)
  iray = 0

  for jRay in 1:nray_max
    if (nr_merge[jRay] < 1)
      continue
    end

    iRay = iRay + 1

    # ! position and width in physical space
    ray = rays[iRay, ix, jy, kz]

    ray.x = 0.5 * (xrmxmg[jray] + xrmnmg[jray])
    ray.dxray = xrmxmg[jray] - xrmnmg[jray]

    ray.y = 0.5 * (yrmxmg[jray] + yrmnmg[jray])
    ray.dyray = yrmxmg[jray] - yrmnmg[jray]

    ray.z = 0.5 * (zrmxmg[jray] + zrmnmg[jray])
    ray.dzray = zrmxmg[jray] - zrmnmg[jray]

    # ! position and width in wavenumber space

    ray.k = 0.5 * (krmxmg[jray] + krmnmg[jray])
    ray.dkray = krmxmg[jray] - krmnmg[jray]

    ray.l = 0.5 * (lrmxmg[jray] + lrmnmg[jray])
    ray.dlray = lrmxmg[jray] - lrmnmg[jray]

    ray.m = 0.5 * (mrmxmg[jray] + mrmnmg[jray])
    ray.dmray = mrmxmg[jray] - mrmnmg[jray]

    # ! intrinsic frequency

    wnrk = ray.k
    wnrl = ray.l
    wnrm = ray.m

    wnrh = sqrt(wnrk^2 + wnrl^2)

    zr = ray.z

    stratification!(zr, 1, NNr)

    omir =
      branchr * sqrt(NNr * wnrh^2 + f_cor_nd^2 * wnrm^2) / sqrt(wnrh^2 + wnrm^2)

    ray.omega = omir

    # ! phase-space volumes

    ray.area_xk = ray.dxray * ray.dkray
    ray.area_yl = ray.dyray * ray.dlray
    ray.area_zm = ray.dzray * ray.dmray

    # ! wave-action density

    if (sizeX > 1)
      fcpspx = ray.area_xk
    else
      fcpspx = 1.0
    end

    if (sizeY > 1)
      fcpspy = ray.area_yl
    else
      fcpspy = 1.0
    end

    fcpspz = ray.area_zm

    if (cons_merge == "wa")
      # ! wave-action density after merging to be determined
      # ! such that the wave action remains the same,
      # ! hence ...

      ray.dens = wadrmg[jray] / (fcpspx * fcpspy * fcpspz)
    elseif (cons_merge == "en")
      # ! wave-action density after merging to be determined
      # ! such that the wave energy remains the same, hence ...

      ray.dens = wadrmg[jray] / (omir * fcpspx * fcpspy * fcpspz)
    else
      # stop 'wrong cons_merge in merge_rayvol'
      exit()
    end
  end
end

function total_ray_volumes()
    return 1
    # nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

  # ! testb
  # ! print*,'before merging nrvloc =',nrvloc
  # ! teste

  # mpi_reduce(nrvloc, nrvtt0, 1, mpi_integer, mpi_sum, root, comm, ierror)
  # mpi_bcast(nrvtt0, 1, mpi_integer, root, comm, ierror)
  return nvrtt0
end

function merge_rayvol(state)
  nvrtt0 = total_ray_volumes()

  f_cor_nd = state.constants.f_coriolis_dim * state.constants.tref

  for kz in 1:nz, jy in 1:ny, ix in 1:nx
    if (nray(ix, jy, kz) <= nray_max)
      continue
    end

    intervals = wave_numbers_interval()

    dwnrk_mg_n, dwnrk_mg_p =
      log_widths(wnrk_max_n, wnrk_min_n, wnrk_max_p, wnrk_min_p, nxray)
    dwnrl_mg_n, dwnrl_mg_p =
      log_widths(wnrl_max_n, wnrl_min_n, wnrl_max_p, wnrl_min_p, nyray)
    dwnrm_mg_n, dwnrm_mg_p =
      log_widths(wnrm_max_n, wnrm_min_n, wnrm_max_p, wnrm_min_p, nzray)

    generate_merged_rv()

    replace_rayvolume()
  end

  nrvtt1 = total_ray_volumes()

  if (master && nrvtt1 < nrvtt0)
    println("after merging nray =", nrvtt1)
  end

  return
end
