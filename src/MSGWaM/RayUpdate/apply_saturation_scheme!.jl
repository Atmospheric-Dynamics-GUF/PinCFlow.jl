function apply_saturation_scheme!(state::State, dt::AbstractFloat)
    (; testcase) = state.namelists.setting
    apply_saturation_scheme!(state, dt, testcase)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractTestCase,
)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractWKBTestCase,
)
    (; wkb_mode) = state.namelists.wkb
    apply_saturation_scheme!(state, dt, wkb_mode)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::SteadyState,
)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::Union{SingleColumn, MultiColumn},
)
    # TODO first and fourth loop are basically identical.

    (; domain, grid) = state
    (; nray, rays, diffusion) = state.wkb
    (; sizex, sizey) = state.namelists.domain
    (; alpha_sat) = state.namelists.wkb
    (; io, jo, i0, i1, j0, j1, k0, k1, nxx, nyy, nzz) = state.domain
    (; lx, ly, dx, dy, dz, ztfc, jac) = state.grid
    (; rhostrattfc) = state.atmosphere

    mb2 = zeros(nxx, nyy, nzz)
    mb2k2 = zeros(nxx, nyy, nzz)

    # Calculate the integrals.
    for kzrv in k0:k1, jyrv in j0:j1, ixrv in i0:i1
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

            omir = rays.omega[iray, ixrv, jyrv, kzrv]

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

            mb2[ix, jy, kz] +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral1

            integral2 = wnrhs * wnrm^2 / omir * facpsp

            mb2k2[ix, jy, kz] +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral2
        end
    end

    # Calculate the turbulent eddy diffusivity.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if mb2k2[ix, jy, kz] == 0 || mb2[ix, jy, kz] < alpha_sat^2 * n2r^2
            diffusion[ix, jy, kz] = 0
        else
            diffusion[ix, jy, kz] =
                (mb2[ix, jy, kz] - alpha_sat^2 * n2r^2) /
                (2 * dt * mb2k2[ix, jy, kz])
        end
    end

    # Reduce the wave-action density.
    for kzrv in k0:k1, jyrv in j0:j1, ixrv in i0:i1
        for iray in 1:nray[ixrv, jyrv, kzrv]

            # Skip ray volumes with zero wave-action density.
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0.0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

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

            wnrk = rays.k[iray, ixrv, jyrv, kzrv]
            wnrl = rays.l[iray, ixrv, jyrv, kzrv]
            wnrm = rays.m[iray, ixrv, jyrv, kzrv]

            kappa = diffusion[ix, jy, kz]

            rays.dens[iray, ixrv, jyrv, kzrv] *=
                max(0, 1 - dt * 2 * kappa * (wnrk^2 + wnrl^2 + wnrm^2))
        end
    end

    # Compute integrals again for diagnostics (exact repetition of the first
    # loop).
    mb2 .= 0
    for kzrv in k0:k1, jyrv in j0:j1, ixrv in i0:i1
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

            omir = rays.omega[iray, ixrv, jyrv, kzrv]

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

            mb2[ix, jy, kz] +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral1
        end
    end

    # Check if saturation is violated.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if mb2[ix, jy, kz] - alpha_sat^2 * n2r^2 > 1.0E-3 * alpha_sat^2 * n2r^2
            println(
                "Saturation violated at (ix, jy, kz) = (",
                ix,
                ", ",
                jy,
                ", ",
                kz,
                ")",
            )
            println("mb2[ix, jy, kz] = ", mb2[ix, jy, kz])
            println("alpha_sat^2 * n2r^2 = ", alpha_sat^2 * n2r^2)
        end
    end

    # Rmove rays with zero wave-action density.
    remove_rays!(state)

    return
end
