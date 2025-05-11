function propagate_rays!(state::State, dt::AbstractFloat, rkstage::Integer)
    (; testcase) = state.namelists.setting
    propagate_rays!(state, dt, rkstage, testcase)
    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    testcase::AbstractTestCase,
)
    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    testcase::AbstractWKBTestCase,
)
    (; wkb_mode) = state.namelists.wkb
    propagate_rays!(state, dt, rkstage, wkb_mode)
    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::AbstractWKBMode,
)
    (; testcase) = state.namelists.setting
    (; branchr, zmin_wkb_dim) = state.namelists.wkb
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; lref, tref) = state.constants
    (; nray_max, nray, cgx_max, cgy_max, cgz_max, rays) = state.wkb
    (; dxray, dyray, dzray, dkray, dlray, dmray, ddxray, ddyray, ddzray) =
        state.wkb.increments
    (; alphark, betark, stepfrac) = state.time
    (; lz, ztildetfc) = state.grid
    (; ko, k0, k1, j0, j1, i0, i1) = state.domain

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    kz0 = testcase == WKBMountainWave() && ko == 0 ? k0 - 1 : k0
    kz1 = k1

    # Initialize RK tendencies at the first RK stage.
    if rkstage == 1
        dxray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        dyray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        dzray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        dkray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        dlray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        dmray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        ddxray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        ddyray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
        ddzray[1:nray_max, i0:i1, j0:j1, kz0:kz1] .= 0.0
    end

    cgx_max[] = 0.0
    cgy_max[] = 0.0
    cgz_max[i0:i1, j0:j1, kz0:kz1] .= 0.0

    for kz in kz0:kz1, jy in j0:j1, ix in i0:i1
        nskip = 0
        for iray in 1:nray[ix, jy, kz]
            (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            (dxr, dyr, dzr) = get_physical_extent(rays, (iray, ix, jy, kz))
            (axk, ayl, azm) = get_surfaces(rays, (iray, ix, jy, kz))

            xr1 = xr - dxr / 2
            xr2 = xr + dxr / 2
            yr1 = yr - dyr / 2
            yr2 = yr + dyr / 2
            zr1 = zr - dzr / 2
            zr2 = zr + dzr / 2

            khr = sqrt(kr^2 + lr^2)

            # Skip ray volumes that have left the domain.
            if testcase != WKBMountainWave()
                if zr1 < ztildetfc[ix, jy, k0 - 2]
                    nskip += 1
                    continue
                end
            end

            n2r1 = interpolate_stratification(zr1, state, N2())
            n2r = interpolate_stratification(zr, state, N2())
            n2r2 = interpolate_stratification(zr2, state, N2())

            omir1 =
                branchr * sqrt(n2r1 * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            omir2 =
                branchr * sqrt(n2r2 * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            # Compute intrinsic zonal group velocity.
            if (sizex > 1)
                cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            end

            # Compute intrinsic meridional group velocity.
            if (sizey > 1)
                cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            end

            # Compute intrinsic vertical group velocities at the vertical edges.
            cgirz1 = -mr * (omir1^2 - fc^2) / (omir1 * (khr^2 + mr^2))
            cgirz2 = -mr * (omir2^2 - fc^2) / (omir2 * (khr^2 + mr^2))

            #-------------------------------
            #      Change of position
            #-------------------------------

            # Update zonal position.

            if sizex > 1 && kz >= k0 && wkb_mode != SingleColumn()
                uxr1 = interpolate_mean_flow(xr1, yr, zr, state, U())
                uxr2 = interpolate_mean_flow(xr2, yr, zr, state, U())

                cgrx1 = cgirx + uxr1
                cgrx2 = cgirx + uxr2

                cgrx = (cgrx1 + cgrx2) / 2

                f = cgrx
                dxray[iray, ix, jy, kz] =
                    dt * f + alphark[rkstage] * dxray[iray, ix, jy, kz]
                rays.x[iray, ix, jy, kz] +=
                    betark[rkstage] * dxray[iray, ix, jy, kz]

                cgx_max[] = max(cgx_max[], abs(cgrx))
            end

            # Update meridional position.

            if sizey > 1 && kz >= k0 && wkb_mode != SingleColumn()
                vyr1 = interpolate_mean_flow(xr, yr1, zr, state, V())
                vyr2 = interpolate_mean_flow(xr, yr2, zr, state, V())

                cgry1 = cgiry + vyr1
                cgry2 = cgiry + vyr2

                cgry = (cgry1 + cgry2) / 2

                f = cgry
                dyray[iray, ix, jy, kz] =
                    dt * f + alphark[rkstage] * dyray[iray, ix, jy, kz]
                rays.y[iray, ix, jy, kz] +=
                    betark[rkstage] * dyray[iray, ix, jy, kz]

                cgy_max[] = max(cgy_max[], abs(cgry))
            end

            # Update vertical position.

            cgrz1 = cgirz1
            cgrz2 = cgirz2

            cgrz = (cgrz1 + cgrz2) / 2

            f = cgrz
            dzray[iray, ix, jy, kz] =
                dt * f + alphark[rkstage] * dzray[iray, ix, jy, kz]
            rays.z[iray, ix, jy, kz] +=
                betark[rkstage] * dzray[iray, ix, jy, kz]

            cgz_max[ix, jy, kz] = max(cgz_max[ix, jy, kz], abs(cgrz))

            # Refraction is only allowed above lz[1] + zmin_wkb_dim / lref.

            if zr > lz[1] + zmin_wkb_dim / lref

                #-------------------------------
                #      Change of wavenumber
                #-------------------------------

                dudxr = interpolate_mean_flow(xr, yr, zr, state, DUDX())
                dudyr = interpolate_mean_flow(xr, yr, zr, state, DUDY())
                dudzr = interpolate_mean_flow(xr, yr, zr, state, DUDZ())

                dvdxr = interpolate_mean_flow(xr, yr, zr, state, DVDX())
                dvdyr = interpolate_mean_flow(xr, yr, zr, state, DVDY())
                dvdzr = interpolate_mean_flow(xr, yr, zr, state, DVDZ())

                dn2dzr = interpolate_stratification(zr, state, DN2DZ())

                dkdt = -dudxr * kr - dvdxr * lr
                dldt = -dudyr * kr - dvdyr * lr
                dmdt =
                    -dudzr * kr - dvdzr * lr -
                    khr^2 * dn2dzr / (2 * omir + (khr^2 + mr^2))

                dkray[iray, ix, jy, kz] =
                    dt * dkdt + alphark[rkstage] * dkray[iray, ix, jy, kz]
                dlray[iray, ix, jy, kz] =
                    dt * dldt + alphark[rkstage] * dlray[iray, ix, jy, kz]
                dmray[iray, ix, jy, kz] =
                    dt * dmdt + alphark[rkstage] * dmray[iray, ix, jy, kz]

                rays.k[iray, ix, jy, kz] +=
                    betark[rkstage] * dkray[iray, ix, jy, kz]
                rays.l[iray, ix, jy, kz] +=
                    betark[rkstage] * dlray[iray, ix, jy, kz]
                rays.m[iray, ix, jy, kz] +=
                    betark[rkstage] * dmray[iray, ix, jy, kz]

                #-------------------------------
                #      Change of extents
                #-------------------------------

                # Update extents in x and k.

                if (sizex > 1 && kz >= k0 && wkb_mode != SingleColumn())
                    ddxdt = cgrx2 - cgrx1

                    ddxray[iray, ix, jy, kz] =
                        dt * ddxdt + alphark[rkstage] * ddxray[iray, ix, jy, kz]

                    rays.dxray[iray, ix, jy, kz] +=
                        betark[rkstage] * ddxray[iray, ix, jy, kz]

                    if rays.dxray[iray, ix, jy, kz] <= 0
                        rays.dxray[iray, ix, jy, kz] *= -1
                    end

                    rays.dkray[iray, ix, jy, kz] =
                        axk / rays.dxray[iray, ix, jy, kz]
                end

                # Update extents in y and l.

                if (sizey > 1 && kz >= k0 && wkb_mode != SingleColumn())
                    ddydt = cgry2 - cgry1

                    ddyray[iray, ix, jy, kz] =
                        dt * ddydt + alphark[rkstage] * ddyray[iray, ix, jy, kz]

                    rays.dyray[iray, ix, jy, kz] +=
                        betark[rkstage] * ddyray[iray, ix, jy, kz]

                    if rays.dyray[iray, ix, jy, kz] <= 0
                        rays.dyray[iray, ix, jy, kz] *= -1
                    end

                    rays.dlray[iray, ix, jy, kz] =
                        ayl / rays.dyray[iray, ix, jy, kz]
                end

                # Update extents in z and m.

                ddzdt = cgrz2 - cgrz1

                ddzray[iray, ix, jy, kz] =
                    dt * ddzdt + alphark[rkstage] * ddzray[iray, ix, jy, kz]

                rays.dzray[iray, ix, jy, kz] +=
                    betark[rkstage] * ddzray[iray, ix, jy, kz]

                if rays.dzray[iray, ix, jy, kz] <= 0
                    rays.dzray[iray, ix, jy, kz] *= -1
                end

                rays.dmray[iray, ix, jy, kz] =
                    azm / rays.dzray[iray, ix, jy, kz]
            end
        end

        if nskip > 0
            println(
                nskip,
                " out of ",
                nray[ix, jy, kz],
                " ray volumes have been skipped in propagate_rays!!",
            )
            println("")
        end
    end

    #-------------------------------
    #     Change of wave action
    #-------------------------------

    if spongelayer && unifiedsponge
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            for iray in 1:nray[ix, jy, kz]
                (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
                alphasponge = 2 * interpolate_sponge(xr, yr, zr, state)
                betasponge = 1 / (1 + alphasponge * stepfrac[rkstage] * dt)
                rays.dens[iray, ix, jy, kz] *= betasponge
            end
        end
    end

    if testcase == WKBMountainWave()
        activate_orographic_source!(state, stepfrac[rkstage] * dt)
    end

    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::SteadyState,
)
    (; sizex, sizey) = state.namelists.domain
    (; testcase) = state.namelists.setting
    (; coriolis_frequency) = state.namelists.atmosphere
    (; spongelayer, unifiedsponge) = state.namelists.sponge
    (; branchr, lsaturation, alpha_sat) = state.namelists.wkb
    (; stepfrac) = state.time
    (; tref) = state.constants
    (; comm, sizezz, nzz, nx, ny, ko, k0, k1, j0, j1, i0, i1, down, up) =
        state.domain
    (; dx, dy, dz, ztildetfc, ztfc, jac) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; u, v) = state.variables.predictands
    (; nray, rays) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    if testcase == WKBMountainWave()
        activate_orographic_source!(state, stepfrac[rkstage] * dt)
    end

    if ko != 0
        nray_down = zeros(Int, nx, ny)
        MPI.Recv!(nray_down, comm; source = down)
        nray[i0:i1, j0:j1, k0 - 1] .= nray_down

        count = maximum(nray[i0:i1, j0:j1, k0 - 1])
        if count > 0
            fields = fieldnames(Rays)
            rays_down = zeros(length(fields), count, nx, ny)
            MPI.Recv!(rays_down, comm; source = down)
            for (index, field) in enumerate(fields)
                @views getfield(rays, field)[1:count, i0:i1, j0:j1, k0 - 1] .=
                    rays_down[index, :, :, :]
            end
        end
    end

    # Loop over grid cells.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1

        # Set the ray-volume count.
        nray[ix, jy, kz] = nray[ix, jy, kz - 1]

        # Set up the saturation integrals.
        integral1 = 0.0
        integral2 = 0.0
        m2b2 = 0.0
        m2b2k2 = 0.0

        # Loop over ray volumes.
        for iray in 1:nray[ix, jy, kz]

            # Prepare the ray volume.
            copy_rays!(rays, (iray, ix, jy, kz - 1), (iray, ix, jy, kz))

            # Skip modes with zero wave-action density.
            if rays.dens[iray, ix, jy, kz - 1] == 0
                continue
            end

            # Set the vertical position (and extent).
            rays.z[iray, ix, jy, kz] =
                ztildetfc[ix, jy, kz - 1] +
                (rays.z[iray, ix, jy, kz - 1] - ztildetfc[ix, jy, kz - 2]) /
                jac[ix, jy, kz - 1] * jac[ix, jy, kz]
            rays.dzray[iray, ix, jy, kz] =
                rays.dzray[iray, ix, jy, kz - 1] * jac[ix, jy, kz] /
                jac[ix, jy, kz - 1]

            # Get the horizontal wavenumbers.
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            khr = sqrt(kr^2 + lr^2)

            # Set the reference level.
            kz0 = max(k0, kz - 1)

            # Compute the vertical group velocity at the level below.
            n2r = interpolate_stratification(
                rays.z[iray, ix, jy, kz0],
                state,
                N2(),
            )
            omir =
                -(u[ix, jy, kz0] + u[ix - 1, jy, kz0]) / 2 * kr -
                (v[ix, jy, kz0] + v[ix, jy - 1, kz0]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                mr = rays.m[iray, ix, jy, kz0]
                cgirz0 = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[iray, ix, jy, kz] = 0.0
                continue
            end

            # Compute the local vertical wavenumber and vertical group velocity.
            n2r = interpolate_stratification(
                rays.z[iray, ix, jy, kz],
                state,
                N2(),
            )
            omir =
                -(u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * kr -
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                mr = -branchr * sqrt(khr^2 * (n2r - omir^2) / (omir^2 - fc^2))
                cgirz = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[iray, ix, jy, kz] = 0.0
                continue
            end

            # Set the local vertical wavenumber.
            rays.m[iray, ix, jy, kz] = mr

            # Set the local wave action density.
            if (spongelayer && unifiedsponge)
                (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
                alphasponge = 2 * interpolate_sponge(xr, yr, zr, state)
                rays.dens[iray, ix, jy, kz] =
                    1 / (
                        1 +
                        alphasponge / cgirz *
                        (rays.z[iray, ix, jy, kz] - rays.z[iray, ix, jy, kz0])
                    ) *
                    cgirz0 *
                    rays.dens[iray, ix, jy, kz0] / cgirz
            else
                rays.dens[iray, ix, jy, kz] =
                    cgirz0 * rays.dens[iray, ix, jy, kz0] / cgirz
            end

            # Cycle if the saturation scheme is turned off.
            if !lsaturation
                continue
            end

            # Get the ray volume extents.
            (dxr, dyr, dzr) = get_physical_extent(rays, (iray, ix, jy, kz))
            (dkr, dlr, dmr) = get_spectral_extent(rays, (iray, ix, jy, kz))

            # Compute the phase space factor.
            dzi = min(dzr, jac[ix, jy, kz] * dz)
            facpsp = dzi / jac[ix, jy, kz] / dz * dmr
            if sizex > 1
                dxi = min(dxr, dx)
                facpsp *= dxi / dx * dkr
            end
            if sizey > 1
                dyi = min(dyr, dy)
                facpsp *= dyi / dy * dlr
            end

            # Compute the saturation integrals.
            integral1 = khr^2 * mr^2 / ((khr^2 + mr^2) * omir) * facpsp
            m2b2 +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] *
                integral1 *
                rays.dens[iray, ix, jy, kz]
            integral2 = khr^2 * mr^2 / omir * facpsp
            m2b2k2 +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] *
                integral2 *
                rays.dens[iray, ix, jy, kz] *
                jac[ix, jy, kz] *
                dz / cgirz
        end

        # Compute the diffusion coefficient
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if m2b2k2 == 0 || m2b2 < alpha_sat^2 * n2r^2
            diffusion = 0.0
        else
            diffusion = (m2b2 - alpha_sat^2 * n2r^2) / (2 * m2b2k2)
        end

        # Reduce the wave action density.
        for iray in 1:nray[ix, jy, kz]
            if !lsaturation
                continue
            end
            if rays.dens[iray, ix, jy, kz] == 0
                continue
            end
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            khr = sqrt(kr^2 + lr^2)
            n2r = interpolate_stratification(
                rays.z[iray, ix, jy, kz],
                state,
                N2(),
            )
            omir =
                -(u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * kr -
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                cgirz = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[iray, ix, jy, kz] = 0.0
                continue
            end
            rays.dens[iray, ix, jy, kz] *= max(
                0,
                1 -
                jac[ix, jy, kz] * dz / cgirz * 2 * diffusion * (khr^2 + mr^2),
            )
        end
    end

    if ko + nzz != sizezz
        @views nray_up = nray[i0:i1, j0:j1, k1]
        MPI.Send(nray_up, comm; dest = up)

        count = maximum(nray[i0:i1, j0:j1, k1])
        if count > 0
            fields = fieldnames(Rays)
            rays_up = zeros(length(fields), count, nx, ny)
            for (index, field) in enumerate(fields)
                @views rays_up[index, :, :, :] .=
                    getfield(rays, field)[1:count, i0:i1, j0:j1, k1]
            end
            MPI.Send(rays_up, comm; dest = up)
        end
    end

    return
end
