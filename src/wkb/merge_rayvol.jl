abstract type AbstractWDir end
struct KDir <: AbstractWDir end
struct LDir <: AbstractWDir end
struct MDir <: AbstractWDir end

struct IntervalWidths{T}
    k::IntervalWidth{T}
    l::IntervalWidth{T}
    m::IntervalWidth{T}
end

struct IntervalWidth{T}
    dwn_n::T
    dwn_p::T
IntervalWidth{T}() where {T} = new(T(0), T(0))
end


mutable struct IntervalBound{T}
  max_n::T
  min_n::T
  max_p::T
  min_p::T

  function IntervalBound{T}() where {T}
    return new(zero(T), zero(T), zero(T), zero(T))
  end
end
mutable struct IntervalBounds{T}
  k::Interval{T}
  l::Interval{T}
  m::Interval{T}

  function IntervalBounds{T}() where {T}
    return new(IntervalBound{T}(), IntervalBound{T}(), IntervalBound{T}())
  end
end

function index(rijk, bounds, widths, rays, dir::AbstractWDir, domain)
    wn = wavenumber(rijk, rays, dir)
    area  = warea(rijk, rays, dir)
    fcpsp = area
    ibounds = intervalbounds(dir, intervals)
    iwidths = intervalwidhts(dir, widths)
    nray = nrays(dir, rays)

    if size(domain, dir) > 1
    if (wn < 0.0)
      if (abs(log(-wn / ibounds.max_n) / iwidths.dwn_n) < 1.0)
        ir = nray / 2 - 1
      else
        ir = Int64(log(-wn / ibounds.min_n) / iwidths.dwn_n) + 1
      end
    elseif (wn == 0.0)
      ir = nray / 2
    else
      if (abs(log(wnrk / intervals.k.max_p) / widths.k.dwn_p) < 1.0)
        ir = nray - 1
      else
        ir =
          Int64(log(wn / ibounds.min_p) / iwidths.dwn_p) +
          nray / 2 + 1
      end
    end

    if (ir < 1)
      # print *, 'ERROR in merge_rayvol: ir_k =', ir_k, '< 1'
      # print *, 'wnrk = ', wnrk
      # print *, 'wnrk_min_n = ', wnrk_min_n
      # print *, 'wnrk_max_n = ', wnrk_max_n
      # print *, 'wnrk_min_p = ', wnrk_min_p
      # print *, 'wnrk_max_p = ', wnrk_max_p
      # stop
      exit()
    elseif (ir > nray - 1)
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
    fcpsp = 1.0
    ir = 1
  end # sizeX > 1

  return ir, fcpsp

end
function adjust(i::IntervalBound)
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

function subdivide(wavenumber, i::IntervalBound{T})

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


function wn_intervals(ijk, nray, rays::Rays, nray, domain)
  bounds = IntervalBounds{Float64}()

  for iray in 1:nray[ix, jy, kz]
    rijk = CartesianIndex(iray, ijk)

    # would be great if we could dispatch on sizex
    if (domain.sizex > 1)
      bounds.k = subdivide(rays.k[rijk], bounds.k)
    end

    if (domain.sizey > 1)
      bounds.l = subdivide(rays.l[rijk], bounds.l)
    end

    bounds.m = subdivide(rays.m[rijk], bounds.m)
  end

  if (domain.sizex > 1)
    bounds.k = adjust(bounds.k)
  end

  if (domain.sizey > 1)
    bounds.l = adjust(bounds.l)
  end
  bounds.m = adjust(bounds.m)
  return bounds
end

function log_widths(interval::IntervalBound, nray)
  mg_n = log(interval.max_n / interval.min_n) / (nray / 2 - 1)
  mg_p = log(interval.max_p / interval.min_p) / (nray / 2 - 1)
  return IntervalWidth(mg_n, mg_p)
end

function generate_merged_rv(ijk, rays, nray, intervals, widths, state)
  domain = state.domain
  wkb = state.wkb
  (; nxray, nyray, nzray, nr_merge, cons_merge) = wkb
  for iray in 1:nray[ix, jy, kz]
      rijk = CartesianIndex(iray, ijk)
    wnrk = rays.k[rijk]
    dwnrk = rays.dkray[rijk]

    wnrl = rays.l[rijk]
    dwnrl = rays.dlray[rijk]

    wnrm = rays.m[rijk]
    dwnrm = rays.dmray[rijk]

    xr = rays.x[rijk]
    dxr = rays.dxray[rijk]

    yr = rays.y[rijk]
    dyr = rays.dyray[rijk]

    zr = rays.z[rijk]
    dzr = rays.dzray[rijk]

    axk = rays.area_xk[rijk]
    ayl = rays.area_yl[rijk]
    azm = rays.area_zm[rijk]

    wdr = rays.dens[rijk]
    omir = rays.omega[rijk]


    ir_k = index(rijk, bounds, widths, rays, KDir(), domain)
    ir_l = index(rijk, bounds, widths, rays, LDir(), domain)
    ir_m = index(rijk, bounds, widths, rays, MDir(), domain)

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

  f_cor_nd = state.constants.f_coriolis_dim * state.constants.tref

  for jray in 1:nray_max
    if (nr_merge[jray] < 1)
      continue
    end

    iray += 1

    # ! position and width in physical space
    ray = get_ray(iray, ix, jy, kz, rays)

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

    fcpspx = ifelse(domain.sizex > 1, ray.area_xk, 1.0)
    fcpspy = ifelse(domain.sizey > 1, ray.area_yk, 1.0)
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

  (; nx, ny, nz) = state.domain
  (; nray, rays, nray_max, nxray, nyray, nzray) = state.wkb

  for kz in 1:nz, jy in 1:ny, ix in 1:nx
      ijk = CartesianIndex(ix, jy, kz)
    if (nray[ix, jy, kz] <= nray_max)
      continue
    end

    interval_bounds = wn_intervals(ijk, rays, nray, state.domain)

    k_width= log_widths(interval_bounds.k, nxray)
    l_width = log_widths(interval_bounds.l, nyray)
    m_width = log_widths(interval_bounds.m, nzray)

    w = IntervalWidths(
        k_width,
        l_width,
        m_width
    )
    generate_merged_rv(ix, jy, kz, rays, nray, interval_bounds, w, state)

    replace_rayvolume(rays, state)
  end

  nrvtt1 = total_ray_volumes()

  # if (master && nrvtt1 < nrvtt0)
  #   println("after merging nray =", nrvtt1)
  # end

  return
end
