# TODO: well, this needs a better name...
struct TemporaryRayVolume
  xn::Float64
  xx::Float64
  yn::Float64
  yx::Float64
  zn::Float64
  zx::Float64
  kn::Float64
  kx::Float64
  ln::Float64
  lx::Float64
  mn::Float64
  mx::Float64
  wa::Float64
end

function compute_intervals(ijk, rays, nray, domain)
  k_min_p = 0.0
  k_max_p = 0.0
  l_min_p = 0.0
  l_max_p = 0.0
  m_min_p = 0.0
  m_max_p = 0.0

  k_min_n = 0.0
  k_max_n = 0.0
  l_min_n = 0.0
  l_max_n = 0.0
  m_min_n = 0.0
  m_max_n = 0.0

  sizex, sizey, _ = domainsize(domain)
  for iRay in 1:nray[ijk]
    rijk = CartesianIndex(iRay, ijk)
    k, l, m = wavenumbers(rijk, rays)

    if sizex > 1
      k_min_p, k_max_p, k_min_n, k_max_n =
        min_max(k, k_min_p, k_max_p, k_min_n, k_max_n)
    end

    if sizey > 1
      l_min_p, l_max_p, l_min_n, l_max_n =
        min_max(l, l_min_p, l_max_p, l_min_n, l_max_n)
    end

    m_min_p, m_max_p, m_min_n, m_max_n =
      min_max(m, m_min_p, m_max_p, m_min_n, m_max_n)
  end

  if sizex > 0
    k_min_p, k_max_p, k_min_n, k_max_n =
      adjust_bounds(k_min_p, k_max_p, k_min_n, k_max_n)
  end

  if sizey > 0
    l_min_p, l_max_p, l_min_n, l_max_n =
      adjust_bounds(l_min_p, l_max_p, l_min_n, l_max_n)
  end

  m_min_p, m_max_p, m_min_n, m_max_n =
    adjust_bounds(m_min_p, m_max_p, m_min_n, m_max_n)

  if sizex > 1
    dk_mg_n = log(k_max_n / k_min_n) / (nxray / 2 - 1)
    dk_mg_p = log(k_max_p / k_min_p) / (nxray / 2 - 1)
  end

  if sizey > 1
    dl_mg_n = log(l_max_n / l_min_n) / (nyray / 2 - 1)
    dl_mg_p = log(l_max_p / l_min_p) / (nyray / 2 - 1)
  end
  dm_mg_n = log(m_max_n / m_min_n) / (nzray / 2 - 1)
  dm_mg_p = log(m_max_p / m_min_p) / (nzray / 2 - 1)

  interval_k =
    (min_n = k_min_n, max_n = k_max_n, mg_n = dk_mg_n, mg_p = dk_mg_p)
  interval_l =
    (min_n = l_min_n, max_n = l_max_n, mg_n = dl_mg_n, mg_p = dl_mg_p)
  interval_m =
    (min_n = m_min_n, max_n = m_max_n, mg_n = dm_mg_n, mg_p = dm_mg_p)

  return interval_k, interval_l, interval_m
end

function adjust_bounds(min_p, max_p, min_n, max_n)
  if min_n == 0 && max_n == 0
    if min_p != 0 && max_p != 0
      min_n = min_p
      max_n = max_p
    else
      # all limits zero only applies if all wnrm = 0
      # hence, just in order to provide some numbers ...
      min_n = 1.0
      max_n = 2.0
    end
  end

  if min_p == 0 && max_p == 0
    if min_n != 0 && max_n != 0
      min_p = min_n
      max_p = max_n
    else
      min_p = 1.0
      max_p = 2.0
    end
  end

  # in order to prevent zero-width intervals ...
  if min_n == max_n
    min_n = 0.5 * min_n
    max_n = 2.0 * max_n
  end
  if min_p == max_p
    min_p = 0.5 * min_p
    max_p = 2.0 * max_p
  end
  return min_p, max_p, min_n, max_n
end

function min_max(wnr, min_p, max_p, min_n, max_n)
  if wnr > 0
    if min_p == 0
      min_p = wnr
    else
      min_p = min(min_p, wnr)
    end
    max_p = max(max_p, wnr)
  elseif (wnr < 0)
    if min_n == 0
      min_n = -wnr
    else
      min_n = min(min_n, -wnr)
    end
    max_n = max(max_n, -wnr)
  end

  return min_p, max_p, min_n, max_n
end

function merge_rayvol!(rays, state::AbstractState)
  (; sizex, sizey, i0, i1, j0, j1, k0, k1) = state.domain
  (; nray, nray_max, rays) = state.wkb

  @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
  nray_before = MPI.Allreduce(nray_before, +, comm)

  nr_merge = Array{TemporaryRayVolume}(undef, nray_max)

  for kz in k0:k1, jy in j0:j1, ix in i0:i1
    ijk = CartesianIndex(ix, jy, kz)
    if nRay[ijk] <= nray_max
      continue
    end

    interval_k, interval_l, interval_m =
      compute_intervals(ijk, rays, nray, state.domain)

    nr_merge .= 0
    # generate merged rvs
    for iray in 1:nray[ijk]
      rijk = CartesianIndex(iray, ijk)
      k, l, m = wavenumbers(rijk, rays)
      x, y, z = positions(rijk, rays)
      axk, ayl, azm = areas(rijk, rays)

      if sizex > 1
        fcpspx = axk
        ir_k = index_klm(rijk, k, interval_k, nxray)
      else
        fcpspx = 1.0
        ir_k = 1
      end
      if sizey > 1
        fcpspx = ayl
        ir_l = index_klm(rijk, l, interval_l, nyray)
      else
        fcpspy = 1.0
        ir_l = 1
      end

      fcpspz = azm
      ir_m = index_klm(rijk, m, interval_m, nzray)
      jray = ray_index(ir_k, ir_l, ir_m, nxray, nyray, nzray, sizex, sizey)

      nr_merge[jray] += 1
      if nr_merge[jray] == 1
        if cons_merge == "wa"
          wa = wdr * fxpspx * fcpspy * fcpspz
        else
          wa = wdr * omir * fcpspx * fcpspy * fcpspz
        end

        nr_merge[jray] = TemporaryRayVolume(
          xr - dxr / 2.0,
          xr + dxr / 2.0,
          yr - dyr / 2.0,
          yr + dyr / 2.0,
          zr - dzr / 2.0,
          zr + dzr / 2.0,
          k - dk / 2.0,
          k + dk / 2.0,
          l - dl / 2.0,
          l + dl / 2.0,
          m - dm / 2.0,
          m + dm / 2.0,
          wa,
        )
      else
        vol = nr_merge[jray]
        vol.xn = min(vol.xn, x - dx / 2.0)
        vol.xx = min(vol.xx, x + dx / 2.0)
        vol.yn = min(vol.yn, y - dy / 2.0)
        vol.yx = min(vol.yx, y + dy / 2.0)
        vol.zn = min(vol.zn, z - dz / 2.0)
        vol.zx = min(vol.zx, z + dz / 2.0)
        vol.kn = min(vol.kn, k - dk / 2.0)
        vol.kx = min(vol.kx, k + dk / 2.0)
        vol.ln = min(vol.ln, l - dl / 2.0)
        vol.lx = min(vol.lx, l + dl / 2.0)
        vol.mn = min(vol.mn, m - dm / 2.0)
        vol.mx = min(vol.mx, m + dm / 2.0)
        if cons_merge == "wa"
          vol.wa += wdr * fcpspx * fcpspy * fcpspz
        elseif cons_merge == "en"
          vol.wa += wdr * omir * fcpspx * fcpspy * fcpspz
        end
      end
    end

    # now replace old ray volumes by merged ones
    #
    iray = 0
    for jray in 1:nray_max
      if nr_merge[jray] < 1
        continue
      end
      iray += 1
      rijk = CartesianIndex(iray, ijk)
      tvol = nr_merge[jray]
      rays.x[rijk] = (tvol.xn + tvol.xx) / 2.0
      rays.y[rijk] = (tvol.yn + tvol.yx) / 2.0
      rays.z[rijk] = (tvol.zn + tvol.zx) / 2.0

      rays.k[rijk] = (tvol.kx + tvol.kx) / 2.0
      rays.l[rijk] = (tvol.ln + tvol.lx) / 2.0
      rays.m[rijk] = (tvol.mn + tvol.mx) / 2.0

      rays.dxray[rijk] = tvol.xx - tvol.xn
      rays.dyray[rijk] = tvol.yy - tvol.yn
      rays.dzray[rijk] = tvol.zz - tvol.zn
      rays.dkray[rijk] = tvol.kx - tvol.kn
      rays.dlray[rijk] = tvol.lx - tvol.ln
      rays.dmray[rijk] = tvol.mx - tvol.mn

      rays.area_xk[rijk] = rays.dxray[rijk] * rays.dkray[rijk]
      rays.area_yl[rijk] = rays.dyray[rijk] * rays.dlray[rijk]
      rays.area_zm[rijk] = rays.dzray[rijk] * rays.dmray[rijk]
      omir = intrinsic_frequency(rijk, branchr, rays)
      rays.omega[rijk] = omir
      # wave action density
      # TODO: as functions
      fcpspx = ifelse(sizex > 1, rays.area_xk[rijk], 1.0)
      fcpspy = ifelse(sizey > 1, rays.area_yl[rijk], 1.0)
      fcpspz = rays.area_zm[rijk]

      if cons_merge == "wa"
        wa = tvol.wa / (fcpspx * fcpspy * fcpspz)
      else
        wa = tvol.wa / (omir * fcpspx * fcpspy * fcpspz)
      end
      rays.dens[rijk] = wa
    end
  end

  # total number of rays after merge
  @views nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
  nray_after = MPI.Allreduce(nray_before, +, comm)
  return nothing
end

function intrinsic_frequency(rijk, branchr, rays)
  k, l, m = wavenumbers(rijk, rays)
  h = sqrt(k^2 + l^2)
  z = rays.z[rijk]
  NNr = stratification(z, 1)

  omir = branchr * sqrt(NNr * h^2 + f_cor_nd^2 * m^2) / sqrt(h^2 + m^2)
  return omir
end

function ray_index(ir_k, ir_l, ir_m, nxray, nyray, nzray, sizex, sizey)
  if sizex > 1
    if sizey > 1
      jray =
        (ir_m - 1) * (nyray - 1) * (nxray - 1) + (ir_l - 1) * (nxray - 1) + ir_k
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
  return jray
end

function index_klm(rijk, w, interval, n_ray)
  if w < 0
    if abs(log(-w / interval.max_n) / interval.mg_n) < 1.0
      ir = n_ray / 2 - 1
    else
      ir = Int64(log(-w / interval.min_n) / interval.mg_n) + 1
    end
  elseif w == 0
    ir = n_ray / 2
  else
    if abs(log(w / interval.max_p) / interval.mg_p < 1.0)
      ir = n_ray - 1
    else
      ir = Int64(log(w / interval.min_p) / interval.mg_p) + n_ray / 2 + 1
    end
  end
  return ir
end
