function wavenumbers(rijk::CartesianIndex{4}, rays::Rays)
    return rays.k[rijk], rays.l[rijk], rays.m[rijk]
end

function positions(rijk::CartesianIndex{4}, rays::Rays)
    return rays.x[rijk], rays.y[rijk], rays.z[rijk]
end

const Bound = SVector{2, Float64}

function update_bound(r, dr, bound)
    return Bound(bound[1] - dr / 2.0, bound[2] + dr / 2.0)
end

function _new_bound(r, dr)
    return Bound(r - dr / 2.0, r + dr / 2.0)
end

struct BinnedRayVolumeBounds
    x::Bound
    y::Bound
    z::Bound
    k::Bound
    l::Bound
    m::Bound
    wa::Float64

    function BinnedRayVolumeBounds(x, dx, y, dy, z, dz, k, dk, l, dl, m, dm, wa)
        bounds = new(
            _new_bound(x, dx),
            _new_bound(y, dy),
            _new_bound(z, dz),
            _new_bound(k, dk),
            _new_bound(l, dl),
            _new_bound(m, dm),
            wa,
        )
        return bounds
    end
end

function midpoint(b::Bound)
    return (b[1] + b[2]) / 2.0
end

function diff(b::Bound)
    return b[1] - b[2]
end

function generate_merged_rvs(nr_merge, ijk, intervals, rays, domain)
    sizex, sizey, _ = domainsize(domain)
    nr_merge .= 0
    # generate merged rvs
    for iray in 1:nray[ijk]
        rijk = CartesianIndex(iray, ijk)
        k, l, m = wavenumbers(rijk, rays)
        x, y, z = positions(rijk, rays)
        dx, dy, dz = position_extents(rijk, rays)
        dk, dl, dm = wavenumber_extents(rijk, rays)
        axk, ayl, azm = areas(rijk, rays)

        if sizex > 1
            fcpspx = axk
            ir_k = index_klm(rijk, k, intervals[1], nxray)
        else
            fcpspx = 1.0
            ir_k = 1
        end
        if sizey > 1
            fcpspx = ayl
            ir_l = index_klm(rijk, l, intervals[2], nyray)
        else
            fcpspy = 1.0
            ir_l = 1
        end

        fcpspz = azm
        ir_m = index_klm(rijk, m, intervals[3], nzray)
        jray = ray_index(ir_k, ir_l, ir_m, nxray, nyray, nzray, sizex, sizey)

        nr_merge[jray] += 1
        if nr_merge[jray] == 1
            if merge_mode == ConstantWaveAction()
                wa = wdr * fxpspx * fcpspy * fcpspz
            elseif merge_mode == ConstantWaveEnergy()
                wa = wdr * omir * fcpspx * fcpspy * fcpspz
            end

            nr_merge[jray] = BinnedRayVolumeBounds(
                x,
                dx,
                y,
                dy,
                z,
                dz,
                k,
                dk,
                l,
                dl,
                m,
                dm,
                wa,
            )
        else
            vol = nr_merge[jray]
            vol.x = update_bound(x, dx, vol.x)
            vol.y = update_bound(y, dy, vol.y)
            vol.z = update_bound(z, dz, vol.z)
            vol.l = update_bound(l, dl, vol.l)
            vol.k = update_bound(k, dk, vol.k)
            vol.m = update_bound(m, dm, vol.m)
            nr_merge[jray] = vol # not sure if this is needed
            if merge_mode == ConstantWaveAction()
                vol.wa += wdr * fcpspx * fcpspy * fcpspz
            elseif merge_mode == ConstantWaveEnergy()
                vol.wa += wdr * omir * fcpspx * fcpspy * fcpspz
            end
        end
    end
end

function replace_rayvolumes!(rays, ijk, nr_merge, nray_max)
    iray = 0
    for jray in 1:nray_max
        if nr_merge[jray] < 1
            continue
        end
        iray += 1
        rijk = CartesianIndex(iray, ijk)
        tvol = nr_merge[jray]
        rays.x[rijk] = midpoint(tvol.x)
        rays.y[rijk] = midpoint(tvol.y)
        rays.z[rijk] = midpoint(tvol.z)

        rays.k[rijk] = midpoint(tvol.k)
        rays.l[rijk] = midpoint(tvol.l)
        rays.m[rijk] = midpoint(tvol.m)

        rays.dxray[rijk] = diff(tvol.x)
        rays.dyray[rijk] = diff(tvol.y)
        rays.dzray[rijk] = diff(tvol.z)

        rays.dkray[rijk] = diff(tvol.k)
        rays.dlray[rijk] = diff(tvol.l)
        rays.dmray[rijk] = diff(tvol.m)

        rays.area_xk[rijk] = rays.dxray[rijk] * rays.dkray[rijk]
        rays.area_yl[rijk] = rays.dyray[rijk] * rays.dlray[rijk]
        rays.area_zm[rijk] = rays.dzray[rijk] * rays.dmray[rijk]
        omir = intrinsic_frequency(rijk, branchr, rays)
        rays.omega[rijk] = omir
        fcpspx = ifelse(sizex > 1, rays.area_xk[rijk], 1.0)
        fcpspy = ifelse(sizey > 1, rays.area_yl[rijk], 1.0)
        fcpspz = rays.area_zm[rijk]

        # wave action density
        # TODO: as functions
        if merge_mode == ConstantWaveAction()
            wa = tvol.wa / (fcpspx * fcpspy * fcpspz)
        elseif merge_mode == ConstantWaveEnergy()
            wa = tvol.wa / (omir * fcpspx * fcpspy * fcpspz)
        end
        rays.dens[rijk] = wa
    end
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

function merge_rayvol!(rays, state::State)
    (; sizex, sizey, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray, nray_max, rays) = state.wkb

    @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    nr_merge = Array{BinnedRayVolumeBounds}(undef, nray_max)

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        ijk = CartesianIndex(ix, jy, kz)
        if nRay[ijk] <= nray_max
            continue
        end

        volume_intervals = compute_intervals(ijk, rays, nray, state.domain)

        generate_merged_rvs!(nr_merge, ijk, volume_intervals, rays, domain)
        replace_rayvolumes!(rays, ijk, nr_merge, nray_max)
    end

    # total number of rays after merge
    @views nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_after = MPI.Allreduce(nray_before, +, comm)
    return nothing
end
