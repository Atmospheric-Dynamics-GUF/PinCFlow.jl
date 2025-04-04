struct widths{T}
  dwnrk_mg_n::T
  dwnrk_mg_p::T
  dwnrl_mg_n::T
  dwnrl_mg_p::T
  dwnrm_mg_n::T
  dwnrm_mg_p::T
end

mutable struct Interval{T}
  max_n::T
  min_n::T
  max_p::T
  min_p::T

  function Interval{T}() where {T}
    return new(zero(T), zero(T), zero(T), zero(T))
  end
end
mutable struct Intervals{T}
  k::Interval{T}
  l::Interval{T}
  m::Interval{T}

  function Intervals{T}() where {T}
    return new(Interval{T}(), Interval{T}(), Interval{T}())
  end
end

function adjust(i::Interval)
  # TODO: this needs a better name
  if (i.min_n == 0.0 && i.max_n == 0.0)
    if (i.min_p != 0.0 && i.max_p != 0.0)
      i.min_n = i.min_p
      i.max_n = i.max_p
    else
      # ! all limits zero only applies if all wnrk = 0
      # ! hence, just in order to provide some numbers ...
      i.min_n = 1.0
      i.max_n = 2.0
    end
  end
  if (i.min_p == 0.0 && i.max_p == 0.0)
    if (i.min_n != 0.0 && i.max_n != 0.0)
      i.min_p = i.min_n
      i.max_p = i.max_n
    else
      # ! all limits zero only applies if all wnrk = 0
      # ! hence, just in order to provide some numbers ...
      i.min_p = 1.0
      i.max_p = 2.0
    end
  end

  # ! in order to prevent zero-width intervals ...

  if (i.min_n == i.max_n)
    i.min_n = 0.5 * i.min_n
    i.max_n = 2.0 * i.max_n
  end

  if (i.min_p == i.max_p)
    i.min_p = 0.5 * i.min_p
    i.max_p = 2.0 * i.max_p
  end
  return i
end

function subdivide(wavenumber, i::Interval)

  # ! subdivide the interval of wavenumbers to be processed
  # ! into up to three ranges:
  # ! (1) if there are negative wavenumbers wnrk < 0, an
  # ! interval -wnrk_max_n <= wnrk <= -wnrk_min_n
  # ! (2) if there are zero wavnumbers wnrk = 0, a slot
  # ! for those
  # ! (3) if there are positive wavenumbers wnrk > 0, an
  # ! interval wnrk_min_p <= wnrk <= wnrk_max_p

  w = wavenumber
  if (w > 0.0)
    if (i.min_p == 0.0)
      i.min_p = w
    else
      i.min_p = min(i.min_p, w)
    end

    i.max_p = max(i.max_p, w)
  elseif (w < 0.0)
    if (i.min_n == 0.0)
      i.min_n = -w
    else
      i.min_n = min(i.min_n, -w)
    end

    i.max_n = max(i.max_n, -w)
  end
  return i
end

function wave_numbers_interval(ix, jy, kz, rays, nray, domain)
  intervals = Intervals{Float64}()

  for iray in 1:nray[ix, jy, kz]
    ray = rays[iray, ix, jy, kz]

    if (domain.sizex > 1)
      intervals.k = subdivide(ray.k, intervals.k)
    end

    if (domain.sizey > 1)
      intervals.l = subdivide(ray.l, intervals.l)
    end

    intervals.m = subdivide(ray.m, intervals.m)
  end

  if (domain.sizex > 1)
    intervals.k = adjust(intervals.k)
  end

  if (domain.sizey > 1)
    intervals.l = adjust(intervals.l)
  end
  intervals.m = adjust(intervals.m)
  return intervals
end

function log_widths(interval::Interval, nray)
  mg_n = log(interval.max_n / interval.min_n) / (nray / 2 - 1)
  mg_p = log(interval.max_p / interval.min_p) / (nray / 2 - 1)
  return (mg_n, mg_p)
end

function generate_merged_rv(ix, jy, kz, rays, nray, intervals, widths, state)
  domain = state.domain
  wkb = state.wkb
  (; nxray, nyray, nzray, nr_merge, cons_merge) = wkb
  for iray in 1:nray[ix, jy, kz]
    ray = rays[iray, ix, jy, kz]
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
    if (domain.sizex > 1)
      fcpspx = axk

      # ! nxRay - 1 intervals for k:
      # ! indices 1 ... nxRay/2 - 1 for negative k
      # ! index nxRay/2 for k = 0
      # ! indices nxRay/2 + 1 ... nxRay - 1 for positive k

      if (wnrk < 0.0)
        if (abs(log(-wnrk / intervals.k.max_n) / widths.dwnrk_mg_n) < 1.0)
          ir_k = nxray / 2 - 1
        else
          ir_k = Int64(log(-wnrk / intervals.k.min_n) / widths.dwnrk_mg_n) + 1
        end
      elseif (wnrk == 0.0)
        ir_k = nxray / 2
      else
        if (abs(log(wnrk / intervals.k.max_p) / widths.dwnrk_mg_p) < 1.0)
          ir_k = nxray - 1
        else
          ir_k =
            Int64(log(wnrk / intervals.k.min_p) / widths.dwnrk_mg_p) +
            nxray / 2 +
            1
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
      elseif (ir_k > nxray - 1)
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

    if (domain.sizey > 1)
      fcpspy = ayl

      if (wnrl < 0.0)
        if (abs(log(-wnrl / intervals.l.max_n) / widths.dwnrl_mg_n) < 1.0)
          ir_l = nyray / 2 - 1
        else
          ir_l = Int64(log(-wnrl / intervals.l.min_n) / widths.dwnrl_mg_n) + 1
        end
      elseif (wnrl == 0.0)
        ir_l = nyray / 2
      else
        if (abs(log(wnrl / intervals.l.max_p) / widths.dwnrl_mg_p) < 1.0)
          ir_l = nyray - 1
        else
          ir_l =
            Int64(log(wnrl / intervals.l.min_p) / widths.dwnrl_mg_p) +
            nyray / 2 +
            1
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
      elseif (ir_l > nyray - 1)
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
      if (abs(log(-wnrm / intervals.m.max_n) / widths.dwnrm_mg_n) < 1.0)
        ir_m = nzray / 2 - 1
      else
        ir_m = Int64(log(-wnrm / intervals.m.min_n) / widths.dwnrm_mg_n) + 1
      end
    elseif (wnrm == 0.0)
      ir_m = nzRay / 2
    else
      if (abs(log(wnrm / intervals.m.max_p) / widths.dwnrm_mg_p) < 1.0)
        ir_m = nzray - 1
      else
        ir_m =
          Int64(log(wnrm / intervals.m.min_p) / widths.dwnrm_mg_p) +
          nzray / 2 +
          1
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
    elseif (ir_m > nzray - 1)
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

    if (domain.sizex > 1)
      if (domain.sizey > 1)
        jray =
          (ir_m - 1) * (nyray - 1) * (nxray - 1) +
          (ir_l - 1) * (nxray - 1) +
          ir_k
      else
        jray = (ir_m - 1) * (nxray - 1) + ir_k
      end
    else
      if (sizey > 1)
        jray = (ir_m - 1) * (nyray - 1) + ir_l
      else
        jray = ir_m
      end
    end

    jray = Int64(jray)
    nr_merge[jray] += 1

    if (nr_merge[jray] == 1)
      wkb.xrmnmg[jray] = xr - 0.5 * dxr
      wkb.xrmxmg[jray] = xr + 0.5 * dxr

      wkb.yrmnmg[jray] = yr - 0.5 * dyr
      wkb.yrmxmg[jray] = yr + 0.5 * dyr

      wkb.zrmnmg[jray] = zr - 0.5 * dzr
      wkb.zrmxmg[jray] = zr + 0.5 * dzr

      wkb.krmnmg[jray] = wnrk - 0.5 * dwnrk
      wkb.krmxmg[jray] = wnrk + 0.5 * dwnrk

      wkb.lrmnmg[jray] = wnrl - 0.5 * dwnrl
      wkb.lrmxmg[jray] = wnrl + 0.5 * dwnrl

      wkb.mrmnmg[jray] = wnrm - 0.5 * dwnrm
      wkb.mrmxmg[jray] = wnrm + 0.5 * dwnrm

      if (cons_merge == "wa")
        # ! wave-action density after merging to be determined
        # ! such that the wave action remains the same,
        # ! hence ...

        wkb.wadrmg[jray] = wdr * fcpspx * fcpspy * fcpspz

      elseif (cons_merge == "en")
        # ! wave-action density after merging to be determined
        # ! such that the wave energy remains the same,
        # ! hence ...

        wkb.wadrmg[jray] = wdr * omir * fcpspx * fcpspy * fcpspz

      else
        # stop 'wrong cons_merge in merge_rayvol'
        exit()
      end
    else
      wkb.xrmnmg[jray] = min(wkb.xrmnmg[jray], xr - 0.5 * dxr)
      wkb.xrmxmg[jray] = max(wkb.xrmxmg[jray], xr + 0.5 * dxr)

      wkb.yrmnmg[jray] = min(wkb.yrmnmg[jray], yr - 0.5 * dyr)
      wkb.yrmxmg[jray] = max(wkb.yrmxmg[jray], yr + 0.5 * dyr)

      wkb.zrmnmg[jray] = min(wkb.zrmnmg[jray], zr - 0.5 * dzr)
      wkb.zrmxmg[jray] = max(wkb.zrmxmg[jray], zr + 0.5 * dzr)

      wkb.krmnmg[jray] = min(wkb.krmnmg[jray], wnrk - 0.5 * dwnrk)
      wkb.krmxmg[jray] = max(wkb.krmxmg[jray], wnrk + 0.5 * dwnrk)

      wkb.lrmnmg[jray] = min(wkb.lrmnmg[jray], wnrl - 0.5 * dwnrl)
      wkb.lrmxmg[jray] = max(wkb.lrmxmg[jray], wnrl + 0.5 * dwnrl)

      wkb.mrmnmg[jray] = min(wkb.mrmnmg[jray], wnrm - 0.5 * dwnrm)
      wkb.mrmxmg[jray] = max(wkb.mrmxmg[jray], wnrm + 0.5 * dwnrm)

      if (cons_merge == "wa")
        # ! wave-action density after merging to be determined
        # ! such that the wave action remains the same,
        # ! hence ...

        wkb.wadrmg[jray] = wkb.wadrmg[jray] + wdr * fcpspx * fcpspy * fcpspz
      elseif (cons_merge == "en")
        # ! wave-action density after merging to be determined
        # ! such that the wave energy remains the same,
        # ! hence ...

        wkb.wadrmg[jRay] += wdr * omir * fcpspx * fcpspy * fcpspz
      else
        # stop 'wrong cons_merge in merge_rayvol'
        exit()
      end
    end
  end
  return
end

function replace_rayvolume(rays, state)
  wkb = state.wkb
  (; nr_merge, nray_max, cons_merge) = wkb
  iray = 0

  for jray in 1:nray_max
    if (nr_merge[jray] < 1)
      continue
    end

    iray += 1

    # ! position and width in physical space
    ray = rays[iray, ix, jy, kz]

    ray.x = 0.5 * (wkb.xrmxmg[jray] + wkb.xrmnmg[jray])
    ray.dxray = wkb.xrmxmg[jray] - wkb.xrmnmg[jray]

    ray.y = 0.5 * (wkb.yrmxmg[jray] + wkb.yrmnmg[jray])
    ray.dyray = wkb.yrmxmg[jray] - wkb.yrmnmg[jray]

    ray.z = 0.5 * (wkb.zrmxmg[jray] + wkb.zrmnmg[jray])
    ray.dzray = wkb.zrmxmg[jray] - wkb.zrmnmg[jray]

    # ! position and width in wavenumber space

    ray.k = 0.5 * (wkb.krmxmg[jray] + wkb.krmnmg[jray])
    ray.dkray = wkb.krmxmg[jray] - wkb.krmnmg[jray]

    ray.l = 0.5 * (wkb.lrmxmg[jray] + wkb.lrmnmg[jray])
    ray.dlray = wkb.lrmxmg[jray] - wkb.lrmnmg[jray]

    ray.m = 0.5 * (wkb.mrmxmg[jray] + wkb.mrmnmg[jray])
    ray.dmray = wkb.mrmxmg[jray] - wkb.mrmnmg[jray]

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

    if (domain.sizeX > 1)
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
  (; nx, ny, nz) = state.domain
  (; nray, rays, nray_max, nxray, nyray, nzray) = state.wkb

  for kz in 1:nz, jy in 1:ny, ix in 1:nx
    if (nray[ix, jy, kz] <= nray_max)
      continue
    end

    intervals = wave_numbers_interval(ix, jy, kz, rays, nray, state.domain)

    dwnrk_mg_n, dwnrk_mg_p = log_widths(intervals.k, nxray)
    dwnrl_mg_n, dwnrl_mg_p = log_widths(intervals.l, nyray)
    dwnrm_mg_n, dwnrm_mg_p = log_widths(intervals.m, nzray)

    w = widths(
      dwnrk_mg_n,
      dwnrk_mg_p,
      dwnrl_mg_n,
      dwnrl_mg_p,
      dwnrm_mg_n,
      dwnrm_mg_p,
    )
    generate_merged_rv(ix, jy, kz, rays, nray, intervals, w, state)

    replace_rayvolume(rays, state)
  end

  nrvtt1 = total_ray_volumes()

  # if (master && nrvtt1 < nrvtt0)
  #   println("after merging nray =", nrvtt1)
  # end

  return
end
