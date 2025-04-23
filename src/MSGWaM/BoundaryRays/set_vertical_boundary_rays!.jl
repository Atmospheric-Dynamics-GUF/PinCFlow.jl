function set_vertical_boundary_rays!(
    state::State,
    zboundaries::SolidWallBoundaries,
)
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, ly, lz, dx, dy, ztildetfc) = state.grid
    (; nray, rays) = state.wkb

    for kz in k0:k1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        nrlc = 0
        for iray in 1:nray[ix, jy, kz]
            zr = rays.z[iray, ix, jy, kz]
            dzr = rays.dzray[iray, ix, jy, kz]
            wnrm = rays.m[iray, ix, jy, kz]

            # Cut ray volumes at the upper boundary.
            if zr - 0.5 * dzr > lz[2]
                continue
            end
            if zr + 0.5 * dzr > lz[2]
                rays.dzray[iray, ix, jy, kz] = lz[2] - zr + 0.5 * dzr
                rays.z[iray, ix, jy, kz] =
                    lz[2] - 0.5 * rays.dzray[iray, ix, jy, kz]
            end

            # Reflect ray volumes at the lower boundary.
            xr = rays.x[iray, ix, jy, kz]
            yr = rays.y[iray, ix, jy, kz]
            ixrv = floor(Int, (xr - lx[1]) / dx) + i0 - io
            jyrv = floor(Int, (yr - ly[1]) / dy) + j0 - jo
            if ztildetfc[ixrv, jyrv, k0 - 1] - zr + 0.5 * dzr > eps()
                rays.z[iray, ix, jy, kz] =
                    2.0 * ztildetfc[ixrv, jyrv, k0 - 1] - zr + dzr
                rays.m[iray, ix, jy, kz] = -wnrm
            end

            nrlc += 1
            if nrlc != iray
                copy_rays!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))
            end
        end
        nray[ix, jy, kz] = nrlc
    end

    return
end

function set_vertical_boundary_rays!(
    state::State,
    zboundaries::PeriodicBoundaries,
)
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, ly, lz, dx, dy, ztildetfc) = state.grid
    (; nray, rays) = state.wkb

    for kz in k0:k1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        if nray[ix, jy, kz] <= 0
            continue
        end
        for iray in 1:nray[ix, jy, kz]
            zr = rays.z[iray, ix, jy, kz]

            xr = rays.x[iray, ix, jy, kz]
            yr = rays.y[iray, ix, jy, kz]
            ixrv = floor(Int, (xr - lx[1]) / dx) + i0 - io
            jyrv = floor(Int, (yr - ly[1]) / dy) + j0 - jo
            zsfc = ztildetfc[ixrv, jyrv, k0 - 1]

            if zr < zsfc
                zr = lz[2] + (zr - zsfc) % (lz[2] - zsfc)
            elseif zr > lz[2]
                zr = zsfc + (zr - lz[2]) % (lz[2] - zsfc)
            end

            rays.z[iray, ix, jy, kz] = zr
        end
    end

    @views nray[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), k0 - 1] .=
        nray[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), k1]
    @views nray[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), k1 + 1] .=
        nray[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), k0]

    for jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        if nray[ix, jy, k0 - 1] > 0
            for iray in 1:nray[ix, jy, k0 - 1]
                copy_rays!(rays, (iray, ix, jy, k1), (iray, ix, jy, k0 - 1))
            end
        end

        if nray[ix, jy, k1 + 1] > 0
            for iray in 1:nray[ix, jy, k1 + 1]
                copy_rays!(rays, (iray, ix, jy, k0), (iray, ix, jy, k1 + 1))
            end
        end

        for kz in (k0 - 1):k0
            if nray[ix, jy, kz] > 0
                for iray in 1:nray[ix, jy, kz]
                    zr = rays.z[iray, ix, jy, kz]
                    zrt = zr - lz[2] + lz[1]

                    xr = rays.x[iray, ix, jy, kz]
                    yr = rays.y[iray, ix, jy, kz]
                    ixrv = floor(Int, (xr - lx[1]) / dx) + i0 - io
                    jyrv = floor(Int, (yr - ly[1]) / dy) + j0 - jo
                    if abs(zrt - ztfc[ixrv, jyrv, kz]) <
                       abs(zr - ztfc[ixrv, jyrv, kz])
                        zr = zrt
                    end

                    rays.z[iray, ix, jy, kz] = zr
                end
            end
        end

        for kz in k1:(k1 + 1)
            if nray[ix, jy, kz] > 0
                for iray in 1:nray[ix, jy, kz]
                    zr = rays.z[iray, ix, jy, kz]
                    zrt = zr + lz[2] - lz[1]

                    xr = rays.x[iray, ix, jy, kz]
                    yr = rays.y[iray, ix, jy, kz]
                    ixrv = floor(Int, (xr - lx[1]) / dx) + i0 - io
                    jyrv = floor(Int, (yr - ly[1]) / dy) + j0 - jo
                    if abs(zrt - ztfc[ixrv, jyrv, kz]) <
                       abs(zr - ztfc[ixrv, jyrv, kz])
                        zr = zrt
                    end

                    rays.z[iray, ix, jy, kz] = zr
                end
            end
        end
    end

    return
end
