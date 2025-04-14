function compute_saturation_integrals(
    state::State,
    indices::NTuple{3, <:Integer},
)
    (; domain, grid) = state
    (; sizex, sizey) = state.namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, dz, jac) = grid
    (; rhostrattfc) = state.atmosphere
    (; nray, rays) = state.wkb

    # Get indices.
    (ixrv, jyrv, kzrv) = indices

    # Initialize Integrals.
    mb2 = 0.0
    mb2k2 = 0.0

    # Loop over ray volumes.
    for iray in 1:nray[ixrv, jyrv, kzrv]

        # Skip ray volumes with zero wave-action density.
        if rays.dens[iray, ixrv, jyrv, kzrv] == 0
            continue
        end

        xr = rays.x[iray, ixrv, jyrv, kzrv]
        yr = rays.y[iray, ixrv, jyrv, kzrv]
        zr = rays.z[iray, ixrv, jyrv, kzrv]

        dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
        dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
        dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

        if sizex > 1
            ix = round(Int, (xr - lx[1] - dx / 2) / dx) + i0 - io
        else
            ix = i0
        end

        if sizey > 1
            jy = round(Int, (yr - ly[1] - dy / 2) / dy) + j0 - jo
        else
            jy = j0
        end

        kz = get_next_half_level(ix, jy, zr, domain, grid)

        n2r = interpolate_stratification(zr, state, N2())

        wnrk = rays.k[iray, ixrv, jyrv, kzrv]
        wnrl = rays.l[iray, ixrv, jyrv, kzrv]
        wnrm = rays.m[iray, ixrv, jyrv, kzrv]

        wnrhs = wnrk^2 + wnrl^2

        dwnrk = rays.dkray[iray, ixrv, jyrv, kzrv]
        dwnrl = rays.dlray[iray, ixrv, jyrv, kzrv]
        dwnrm = rays.dmray[iray, ixrv, jyrv, kzrv]

        omir = compute_intrinsic_frequency(state, (iray, ixrv, jyrv, kzrv))

        densr = rays.dens[iray, ixrv, jyrv, kzrv]

        dzi = min(dzr, jac[ix, jy, kz] * dz)
        facpsp = dzi / jac[ix, jy, kz] / dz * dwnrm

        if sizex > 1
            dxi = min(dxr, dx)
            facpsp = facpsp * dxi / dx * dwnrk
        end

        if sizey > 1
            dyi = min(dyr, dy)
            facpsp = facpsp * dyi / dy * dwnrl
        end

        integral1 = wnrhs * wnrm^2 / ((wnrhs + wnrm^2) * omir) * facpsp

        mb2 += 2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral1

        integral2 = wnrhs * wnrm^2 / omir * facpsp

        mb2k2 += 2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral2
    end

    # Return the results.
    return (mb2, mb2k2)
end