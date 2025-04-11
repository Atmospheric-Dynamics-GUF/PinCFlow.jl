function transport_rayvol(state::State, dt::AbstractFloat, rkstage::Integer)
    (; wkb_mode) = state.namelists.wkb
    transport_rayvol(state, dt, rkstage, wkb_mode)
    return
end

function transport_rayvol(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::Union{SingleColumn, MultiColumn},
)
    (; testcase) = state.namelists.setting
    (; branchr, zmin_wkb) = state.namelists.wkb
    (; sizex, sizey) = state.namelists.domain
    (; nray, cgx_max, cgy_max, cgz_max_tfc, rays, f_cor_nd) = state.wkb
    (; dxray, dkray, ddxray) = state.wkb.increments
    (; alphark, betark, stepfrac) = state.time
    lz = state.grid.lz
    (; k0, k1, j0, j1, i0, i1) = state.domain

    if (rkstage == 1)
        # ! initialize RK-tendencies at first RK stage
        dxray .= 0.0
        dkray .= 0.0
        ddxray .= 0.0
    end

    cgx_max = 0.0
    cgy_max = 0.0
    cgz_max = 0.0
    cgz_max_tfc .= 0.0

    kz0 = ifelse(testcase == WKBMountainWave(), 0, 1)

    for kz in kz0:k1, jy in j0:j1, ix in i0:i1
        ijk = CartesianIndex(ix, jy, kz)
        nskip = 0

        if (nray[ijk] < 1)
            continue
        end

        for iray in 1:nray[ijk]
            rijk = CartesianIndex(iray, ijk)
            wnrk, wnrl, wnrm = get_wavenumbers(rijk, rays)
            wnrh = sqrt(wnrk^2 + wnrl^2)
            xr, yr, zr = get_positions(rijk, rays)
            dxr, dyr, dzr = get_physical_extents(rijk, rays)
            xr1 = xr - 0.5 * dxr
            xr2 = xr + 0.5 * dxr
            yr1 = yr - 0.5 * dyr
            yr2 = yr + 0.5 * dyr
            zr1 = zr - 0.5 * dzr
            zr2 = zr + 0.5 * dzr

            # !skip ray volumes that have left the domain
            if testcase != WKBMountainWave()
                if (zr1 < state.grid.ztildetfc[ix, jy, -1])
                    nskip = nskip + 1
                    continue
                end
            end

            nnr = stratification(zr, state, N2())
            nnr1 = stratification(zr1, state, N2())
            nnr2 = stratification(zr2, state, N2())

            omir1 =
                branchr * sqrt(nnr1 * wnrh^2 + f_cor_nd^2 * wnrm^2) /
                sqrt(wnrh^2 + wnrm^2)

            omir =
                branchr * sqrt(nnr * wnrh^2 + f_cor_nd^2 * wnrm^2) /
                sqrt(wnrh^2 + wnrm^2)

            if (nnr2 <= 0.0)
                # print *, 'NNr2 =', NNr2, '<= 0.0 at'
                # print *, 'zr2 =', zr2, 'from'
                # print *, 'ray(iRay,ix,jy,kz)%z =', rays[rijk]%z
                # print *, 'ray(iRay,ix,jy,kz)%dzray =', rays[rijk]%dzray
                # print *, 'iRay,ix,jy,kz =', rijk
                exit()
            end

            if (wnrh <= 0.0)
                # print *, 'wnrh =', wnrh, '<= 0.0 from'
                # print *, 'wnrk =', wnrk
                # print *, 'wnrl =', wnrl
                # print *, 'iRay,ix,jy,kz =', rijk
                exit()
            end

            omir2 =
                branchr * sqrt(nnr2 * wnrh^2 + f_cor_nd^2 * wnrm^2) /
                sqrt(wnrh^2 + wnrm^2)

            rays.omega[rijk] = omir

            # ! intrinsic group velocities at the respective edges of
            # ! the ray volumes

            if (sizex > 1)
                # ! intrinsic group velocity in x direction not depending
                # ! on x
                cgirx = wnrk * (nnr - omir^2) / (omir * (wnrh^2 + wnrm^2))
            end

            if (sizey > 1)
                # ! intrinsic group velocity in y direction not depending
                # ! on y
                cgiry = wnrl * (nnr - omir^2) / (omir * (wnrh^2 + wnrm^2))
            end

            # ! intrinsic vertical group velocity depending on z
            # ! (via the stratification)
            cgirz1 =
                -wnrm * (omir1^2 - f_cor_nd^2) / (omir1 * (wnrh^2 + wnrm^2))
            cgirz2 =
                -wnrm * (omir2^2 - f_cor_nd^2) / (omir2 * (wnrh^2 + wnrm^2))

            ## until here everything is the same for single_column and transient
            if sizex > 1 && kz >= k0 && mode != SingleColumn()
                uxr1 = meanflow(xr1, yr, zr, state, U())
                uxr2 = meanflow(xr2, yr, zr, state, U())

                cgrx1 = cgirx + uxr1
                cgrx2 = cgirx + uxr2

                cgrx = 0.5 * (cgrx1 + cgrx2)

                F = cgrx
                dxray[1, rijk] = dt * F + alphark(rkstage) * dxray(1, rijk)
                rays.x[rijk] += betark(rkstage) * dxray[1, rijk]

                cgx_max = max(cgx_max, abs(cgrx))
            end

            if sizey > 1 && kz >= k0 && mode != SingleColumn()
                vyr1 = meanflow(xr, yr1, zr, state, V())
                vyr2 = meanflow(xr, yr2, zr, state, V())

                cgry1 = cgiry + vyr1
                cgry2 = cgiry + vyr2

                cgry = 0.5 * (cgry1 + cgry2)
                F = cgry
                dxray[2, rijk] = dt * F + alphark(rkstage) * dxray[2, rijk]
                rays.y[rijk] += betark(rkstage) * dxray[2, rijk]

                cgy_max = max(cgy_max, abs(cgry))
            end

            # !-----------------------------
            # !     vertical displacement
            # !-----------------------------

            # ! RK update

            # ! in line with the asymptotic results the vertcal wind is
            # ! NOT added to the intrinsic vertical group velocity
            # ! should one want to change this, one would also have to
            # ! take the vertical-wind gradient into account in the
            # ! prognostic equations for the wave number
            #
            cgrz1 = cgirz1
            cgrz2 = cgirz2

            cgrz = 0.5 * (cgrz1 + cgrz2)

            F = cgrz
            dxray[3, rijk] = dt * F + alphark[rkstage] * dxray[3, rijk]
            rays.z[rijk] += betark[rkstage] * dxray[3, rijk]

            cgz_max = max(cgz_max, abs(cgrz))

            cgz_max_tfc[ix, jy, kz] = max(cgz_max_tfc[ix, jy, kz], abs(cgrz))

            # !-------------------------------
            # !    change of wavenumber
            # !-------------------------------

            # ! wave refraction only above lz(0) + zmin_wkb
            if (zr > lz[0] + zmin_wkb)
                #   # ! RK procedure

                # TODO: we updated xr etc. what does fortran do here? make a copy or reference?
                dudxr = meanflow(xr, yr, zr, state, DUDX())
                dudyr = meanflow(xr, yr, zr, state, DUDY())
                dudzr = meanflow(xr, yr, zr, state, DUDZ())

                dvdxr = meanflow(xr, yr, zr, state, DVDX())
                dvdyr = meanflow(xr, yr, zr, state, DVDY())
                dvdzr = meanflow(xr, yr, zr, state, DVDZ())

                if (zr < lz[0] - dz)
                    # print *, 'ERROR IN setup_wkb: LOWER EDGE OF RAY  VOLUME', &
                    #     &rijk, 'TOO LOW'
                    exit()
                end

                if (zr < lz[0] - dz)
                    # print *, 'ERROR IN transport_rayvol: RAY VOLUME', iRay, ix, &
                    #     &jy, kz, 'TOO LOW'
                    exit()
                end

                dnndzr = stratification(zr, state, DN2DZ())

                dkdt = -dudxr * wnrk - dvdxr * wnrl
                dldt = -dudyr * wnrk - dvdyr * wnrl
                dmdt =
                    -dudzr * wnrk - dvdzr * wnrl -
                    wnrh^2 * dnndzr / (2.0 * omir + (wnrh^2 + wnrm^2))

                dkray[1, rijk] = dt * dkdt + alphark[rkstage] * dkray[1, rijk]
                dkray[2, rijk] = dt * dldt + alphark[rkstage] * dkray[2, rijk]
                dkray[3, rijk] = dt * dmdt + alphark[rkstage] * dkray[3, rijk]

                rays.k[rijk] += betark[rkstage] * dkray[1, rijk]
                rays.l[rijk] += betark[rkstage] * dkray[2, rijk]
                rays.m[rijk] += betark[rkstage] * dkray[3, rijk]

                # !----------------------------------------------
                # !    change of wave-number width of ray volumes
                # !----------------------------------------------

                # ! dk

                if (sizex > 1 && kz > 0 && mode != SingleColumn())
                    ddxdt = cgrx2 - cgrx1

                    ddxray[1, rijk] =
                        dt * ddxdt + alphark[rkstage] * ddxray[1, rijk]

                    rays.dxray[rijk] += betark[rkstage] * ddxray[1, rijk]

                    if (rays.dxray[rijk] <= 0.0)
                        # print *, 'dxray(', rijk, ') <= 0.0  ==> time &
                        #     &step too large?'
                        rays.dxray[rijk] *= -1
                    end

                    rays.dkray[rijk] = rays.area_xk[rijk] / rays.dxray[rijk]
                end

                # ! dl

                if (sizey > 1 && kz > 0 && mode != SingleColumn())
                    ddydt = cgry2 - cgry1

                    ddxray[2, rijk] =
                        dt * ddydt + alphark[rkstage] * ddxray[2, rijk]

                    rays.dyray[rijk] += betark[rkstage] * ddxray[2, rijk]

                    if (rays.dyray[rijk] <= 0.0)
                        # print *, 'dyray(', rijk, ') <= 0.0  ==> time &
                        #     &step too large?'
                        rays.dyray[rijk] *= -1
                    end

                    rays.dlray[rijk] = rays.area_yl[rijk] / rays.dyray[rijk]
                end

                # !dm

                ddzdt = cgrz2 - cgrz1

                ddxray[3, rijk] =
                    dt * ddzdt + alphark[rkstage] * ddxray[3, rijk]

                rays.dzray[rijk] += betark[rkstage] * ddxray[3, rijk]

                if (rays.dzray[rijk] <= 0.0)
                    # print *, 'dzray(', rijk, ') <= 0.0  ==> time step &
                    #     &too large?'
                    rays.dzray[rijk] *= -1
                end

                rays.dmray[rijk] = rays.area_zm[rijk] / rays.dzray[rijk]
            end

            # !-----------------------------------
            # ! update of the intrinsic frequency
            # !-----------------------------------

            wnrk = rays.k[rijk]
            wnrl = rays.l[rijk]
            wnrm = rays.m[rijk]

            wnrh = sqrt(wnrk^2 + wnrl^2)

            zr = rays.z[rijk]

            nnr = stratification(zr, state, N2())
            omir =
                branchr * sqrt(nnr * wnrh^2 + f_cor_nd^2 * wnrm^2) /
                sqrt(wnrh^2 + wnrm^2)

            rays.omega[rijk] = omir
            if (spongelayer && unifiedsponge)
                xr, yr, zr = get_positions(rijk, rays)
                alphasponge = 2.0 * interpolate_sponge(xr, yr, zr)
                betasponge = 1.0 / (1.0 + alphasponge * stepfrac[rkstage] * dt)
                rays.dens[rijk] *= betasponge
            end
        end # ray loop
        if (nskip > 0)
            # print *, nskip, 'r.v. skipped in transport_rayvol out of', &
            # &nRay(ix, jy, kz)
        end
    end # grid loop

    if (testcase == WKBMountainWave() && wkb_mode != SteadyState())
        orographic_source(var, ray, time, stepfrac(rkstage) * dt)
    end

    return
end

function transport_rayvol(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::SteadyState,
)
    (; testcase) = state.namelists.setting
    (; branchr, zmin_wkb) = state.namelists.wkb
    (; sizex, sizey) = state.namelists.domain
    (; nray, cgz_max_tfc, rays, f_cor_nd) = state.wkb
    (; dxray, dkray, ddxray) = state.wkb.increments
    (; alphark, betark, stepfrac) = state.time
    lz = state.grid.lz
    (; k0, k1, j0, j1, i0, i1) = state.domain
    if testcase == WKBMountainWave()
        orographic_source(var, ray, time, stepFrac(RKStage) * dt)
    end

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        # ! Set ray-volume count.
        nray[ix, jy, kz] = nray[ix, jy, kz - 1]

        # ! Set up saturation computation.
        integral1 = 0.0
        integral2 = 0.0
        m2b2 = 0.0
        m2b2k2 = 0.0
        # ! Loop over ray volumes.
        for iray in 1:nray(ix, jy, kz)
            rijk = CartesianIndex(iray, ix, jy, kz)
            # ! Prepare ray volume.
            copy_rayvolume!(rays, (iray, ix, jy, kz), (iray, ix, jy, kz - 1))

            # ! Skip modes with zero wave-action density.
            if rays.dens[iray, ix, jy, kz - 1] == 0
                continue
            end

            # ! Set vertical position (and extent).
            rays.z[rijk] =
                ztildetfc[ix, jy, kz - 1] + rays.z[iray, ix, jy, kz - 1] -
                ztildetfc[ix, jy, kz - 2] / jac[ix, jy, kz - 1] *
                jac[ix, jy, kz]
            rays.dzray[iray, ix, jy, kz] =
                rays.dzray[iray, ix, jy, kz - 1] * jac[ix, jy, kz] /
                jac[ix, jy, kz - 1]

            # ! Get horizontal wavenumbers.
            k, l, _ = wavenumbers(rijk, rays)
            wnrh = sqrt(k^2.0 + l^2.0)

            # ! Set reference level.
            kz0 = max(1, kz - 1)
            rijk0 = CartesianIndex(iray, ix, jy, kz, kz0)
            # ! Compute vertical group velocity at the level below.
            nn_nd = stratification(rays.z[rijk0], state, N2())
            omir = rays.omega[rijk0].omega

            if (branchr * omir > f_cor_nd && branchr * omir < sqrt(nn_nd))
                wnrm = rays.m[rijk0]
                cgirz0 =
                    wnrm * (f_cor_nd^2 - nn_nd) * wnrh^2 / omir /
                    (wnrh^2 + wnrm^2)^2
            else
                rays.dens[rijk0] = 0.0
                rays.dens[rijk] = 0.0
                continue
            end
            # ! Compute local intrinsic frequency, vertical
            # ! wavenumber and vertical group velocity.
            nn_nd = stratification(rays.z[rijk], state, N2())
            omir =
                -0.5 * (var.u[ix, jy, kz] + var.u[ix - 1, jy, kz]) * wnrk -
                0.5 * (var.v[ix, jy, kz] + var.v[ix, jy - 1, kz]) * wnrl
            if (branchr * omir > f_cor_nd && branchr * omir < sqrt(nn_nd))
                wnrm =
                    -branchr *
                    sqrt(wnrh^2 * (nn_nd - omir^2) / (omir^2 - f_cor_nd^2))
                cgirz =
                    wnrm * (f_cor_nd^2 - NN_nd) * wnrh^2 / omir /
                    (wnrh^2 + wnrm^2)^2
            else
                rays.dens[rijk0] = 0.0
                rays.dens[rijk] = 0.0
                continue
            end

            # ! Set local intrinsic frequency and vertical wavenumber.
            rays.omega[rijk] = omir
            rays.m[rijk] = wnrm

            # ! Set local wave action density.
            if (spongeLayer && unifiedSponge)
                xr = rays.x[rijk]
                yr = rays.y[rijk]
                zr = rays.z[rijk]
                alphasponge = 2.0 * interpolate_sponge(xr, yr, zr)
                rays.dens[rijk] =
                    1.0 / (
                        1.0 +
                        alphaSponge / cgirz * (rays.z[rijk].z - rays.z[rijk0])
                    ) *
                    cgirz0 *
                    rays.dens[iray, ix, jy, kz0] / cgirz
            else
                rays.dens[rijk] = cgirz0 * rays.dens[iray, ix, jy, kz0] / cgirz
            end
            # ! Cycle if saturation scheme is turned off.
            if (!lsaturation)
                continue
            end

            # ! Get ray volume extents.
            dxr, dyr, dzr = get_physical_extents(rijk, rays)
            dwnrk, dwnrl, dwnrm = get_spectral_extents(rijk, rays)

            # ! Compute phase space factor.
            dzi = min(dzr, jac[ix, jy, kz] * dz)
            facpsp = dzi / jac[ix, jy, kz] / dz * dwnrm

            if (sizex > 1)
                dxi = min(dxr, dx)
                facpsp = facpsp * dxi / dx * dwnrk
            end
            if (sizey > 1)
                dyi = min(dyr, dy)
                facpsp = facpsp * dyi / dy * dwnrl
            end

            # ! Update saturation amplitude.
            integral1 = wnrh^2 * wnrm^2 / ((wnrh^2 + wnrm^2) * omir) * facpsp
            m2b2 =
                m2b2 +
                2.0 * nn_nd^2 / rhostrattfc[ix, jy, kz] *
                integral1 *
                rays.dens[rijk]

            integral2 = wnrh^2 * wnrm^2 / omir * facpsp
            m2b2k2 =
                m2b2k2 +
                2.0 * nn_nd^2 / rhostrattfc[ix, jy, kz] *
                integral2 *
                rays.dens[rijk] *
                jac[ix, jy, kz] *
                dz / cgirz
        end
        # ! Compute diffusion coefficient
        nn_nd = stratification(ztfc[ix, jy, kz], state, N2())
        if (m2b2k2 == 0.0 || m2b2 < alpha_sat^2 * nn_nd^2)
            diffusion = 0.0
        else
            diffusion = (m2b2 - alpha_sat^2 * nn_nd^2) / (2.0 * m2b2k2)
        end
        # ! Reduce wave action density.
        for iray in 1:nray(ix, jy, kz)
            if (!lsaturation)
                continue
            end
            if (rays.dens[rijk] == 0.0)
                continue
            end
            wnrk = rays.k[rijk]
            wnrl = rays.l[rijk]
            wnrm = rays.m[rijk]
            wnrh = sqrt(wnrk^2 + wnrl^2)
            nn_nd = stratification(rays.z[rijk], state, N2())
            omir = rays.omega[rijk]
            if (branchr * omir > f_cor_nd && branchr * omir < sqrt(nn_nd))
                cgirz =
                    wnrm * (f_cor_nd^2 - nn_nd) * wnrh^2 / omir /
                    (wnrh^2 + wnrm^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[rijk].dens = 0.0
                continue
            end
            rays.dens[rijk] =
                rays[rijk].dens * max(
                    0.0,
                    1.0 -
                    jac[x, jy, kz] * dz / cgirz *
                    2.0 *
                    diffusion *
                    (wnrh^2 + wnrm^2),
                )
        end
    end
    if (sizex > 1)
        setboundary_rayvol_x(ray)
    end
    if (sizey > 1)
        setboundary_rayvol_y(ray)
    end

    return
end
