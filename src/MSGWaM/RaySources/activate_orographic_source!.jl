function activate_orographic_source!(
    state::State,
    omi_ini::AbstractArray{<:AbstractFloat, 4},
    wnk_ini::AbstractArray{<:AbstractFloat, 4},
    wnl_ini::AbstractArray{<:AbstractFloat, 4},
    wnm_ini::AbstractArray{<:AbstractFloat, 4},
    wad_ini::AbstractArray{<:AbstractFloat, 4},
)

    # Get all necessary fields.
    (; f_coriolis_dim) = state.namelists.atmosphere
    (; branchr, blocking, long_threshold, nwm) = state.namelists.wkb
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, ztildetfc, k_spectrum, l_spectrum, topography_spectrum) =
        state.grid
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; u, v) = state.variables.predictands
    (; zb) = state.wkb

    # Set Coriolis parameter.
    f_cor_nd = f_coriolis_dim * tref

    # Iterate over surface grid cells.
    for jy in j0:j1, ix in i0:i1

        # Average mean wind, reference density and buoyancy frequency. This should
        # be done without a vertical loop.
        uavg = 0.0
        vavg = 0.0
        rhoavg = 0.0
        bvsavg = 0.0
        dzsum = 0.0
        for kz in k0:k1
            uavg +=
                0.5 * (u[ix, jy, kz] + u[ix - 1, jy, kz]) * jac[ix, jy, kz] * dz
            vavg +=
                0.5 * (v[ix, jy, kz] + v[ix, jy - 1, kz]) * jac[ix, jy, kz] * dz
            rhoavg += rhostrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            bvsavg += bvsstrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            dzsum += jac[ix, jy, kz] * dz
            if ztildetfc[ix, jy, kz] >
               ztildetfc[ix, jy, k0 - 1] +
               sum(abs.(topography_spectrum[:, ix, jy]))
                break
            end
        end
        uavg = uavg / dzsum
        vavg = vavg / dzsum
        rhoavg = rhoavg / dzsum
        bvsavg = bvsavg / dzsum

        # Determine the blocked layer.
        if blocking && sum(abs(topography_spectrum[:, ix, jy])) > 0.0
            long =
                sqrt(bvsavg) / sqrt(uavg^2.0 + vavg^2.0) *
                sum(abs.(topography_spectrum[:, ix, jy]))
            ratio = min(1.0, long_threshold / long)
            zb[ix, jy] =
                ztildetfc[ix, jy, 0] +
                sum(abs.(topography_spectrum[:, ix, jy])) * (1.0 - 2.0 * ratio)
        elseif blocking
            ratio = 1.0
            zb[ix, jy] = ztildetfc[ix, jy, 0]
        else
            ratio = 1.0
        end

        # Set launch level.
        kz = k0 - 1

        # Iterate over wave modes.
        for iwm in 1:nwm

            # Compute intrinsic frequency, wavenumbers and wave-action density.
            (omi, wnk, wnl, wnm, wad) = compute_orographic_mode(
                ratio * topography_spectrum[iwm, ix, jy],
                k_spectrum[iwm, ix, jy],
                l_spectrum[iwm, ix, jy],
                uavg,
                vavg,
                rhoavg,
                bvsavg,
                f_cor_nd,
                branchr,
            )

            # Save the results.
            omi_ini[iwm, ix, jy, kz] = omi
            wnk_ini[iwm, ix, jy, kz] = wnk
            wnl_ini[iwm, ix, jy, kz] = wnl
            wnm_ini[iwm, ix, jy, kz] = wnm
            wad_ini[iwm, ix, jy, kz] = wad
        end
    end
    return
end

function activate_orographic_source!(state::State, dt::AbstractFloat)
    # Get all necessary fields.
    (; f_cor_nd) = state.namelists.atmosphere
    (;
        nrxl,
        nryl,
        nrzl,
        nrk_init,
        nrl_init,
        nrm_init,
        fac_dk_init,
        fac_dl_init,
        fac_dm_init,
        branchr,
        blocking,
        long_threshold,
        launch_algorithm,
    ) = state.namelists.wkb
    (; tref) = state.constants
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (;
        dx,
        dy,
        dz,
        x,
        y,
        ztfc,
        jac,
        ztildetfc,
        k_spectrum,
        l_spectrum,
        topography_spectrum,
    ) = state.grid
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; u, v) = state.variables.predictands
    (; ir_sfc, ix2_sfc, jy2_sfc, kz2_sfc, ik_sfc, jl_sfc, km_sfc, iwm_sfc) =
        state.wkb.surface_indices
    (; nray_wrk, nray, rays) = state.wkb

    # Set Coriolis parameter.
    f_cor_nd = f_coriolis_dim * tref

    # Iterate over surface grid cells.
    for jy in j0:j1, ix in i0:i1

        # Average mean wind, reference density and buoyancy frequency. This should
        # be done without a vertical loop.
        uavg = 0.0
        vavg = 0.0
        rhoavg = 0.0
        bvsavg = 0.0
        dzsum = 0.0
        for kz in k0:k1
            uavg +=
                0.5 * (u[ix, jy, kz] + u[ix - 1, jy, kz]) * jac[ix, jy, kz] * dz
            vavg +=
                0.5 * (v[ix, jy, kz] + v[ix, jy - 1, kz]) * jac[ix, jy, kz] * dz
            rhoavg += rhostrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            bvsavg += bvsstrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            dzsum += jac[ix, jy, kz] * dz
            if ztildetfc[ix, jy, kz] >
               ztildetfc[ix, jy, k0 - 1] +
               sum(abs.(topography_spectrum[:, ix, jy]))
                break
            end
        end
        uavg = uavg / dzsum
        vavg = vavg / dzsum
        rhoavg = rhoavg / dzsum
        bvsavg = bvsavg / dzsum

        # Determine the blocked layer.
        if blocking && sum(abs(topography_spectrum[:, ix, jy])) > 0.0
            long =
                sqrt(bvsavg) / sqrt(uavg^2.0 + vavg^2.0) *
                sum(abs.(topography_spectrum[:, ix, jy]))
            ratio = min(1.0, long_threshold / long)
            zb[ix, jy] =
                ztildetfc[ix, jy, 0] +
                sum(abs.(topography_spectrum[:, ix, jy])) * (1.0 - 2.0 * ratio)
        elseif blocking
            ratio = 1.0
            zb[ix, jy] = ztildetfc[ix, jy, 0]
        else
            ratio = 1.0
        end

        # Set launch level.
        kz = k0 - 1

        # Loop over surface ray volumes.
        for i_sfc in 1:n_sfc
            iray = ir_sfc[i_sfc, ix, jy]

            # Set surface indices.
            ix2 = ix2_sfc[i_sfc]
            jy2 = jy2_sfc[i_sfc]
            kz2 = kz2_sfc[i_sfc]
            ik = ik_sfc[i_sfc]
            jl = jl_sfc[i_sfc]
            km = km_sfc[i_sfc]
            iwm = iwm_sfc[i_sfc]

            # Compute intrinsic frequency, wavenumbers and wave-action density.
            (omir, wnrk, wnrl, wnrm, wadr) = compute_orographic_mode(
                ratio * topography_spectrum[iwm, ix, jy],
                k_spectrum[iwm, ix, jy],
                l_spectrum[iwm, ix, jy],
                uavg,
                vavg,
                rhoavg,
                bvsavg,
                f_cor_nd,
                branchr,
            )

            # Get vertical position and extent of old ray volume.
            if iray > 0
                zr = rays.z[iray, ix, jy, kz]
                dzr = rays.dzray[iray, ix, jy, kz]
            else
                zr = 0.0
                dzr = 0.0
            end

            # Check if a new ray volume is to be launched and clip the old one.
            # Three cases are distinguished.
            # (1) There is no ray volume with nonzero wave-action density. A new
            #     ray volume is launched.
            # (2) There is a ray volume with nonzero wave-action density, which
            #     has partially passed the lower boundary. It is clipped and the
            #     part below the lower boundary discarded before a new ray
            #     volume is launched.
            # (3) There is a ray volume with nonzero wave-action density, which
            #     has not yet crossed the lower boundary. It is replaced with a
            #     new one.
            if wkb_mode == SteadyState()
                if iray < 0
                    nray[ix, jy, kz] += 1
                    iray = nray[ix, jy, kz]
                    ir_sfc[i_sfc, ix, jy] = iray
                end

                if wadr == 0.0
                    rays.dens[iray, ix, jy, kz] = 0.0
                    continue
                end
            else
                if wadr != 0.0
                    if launch_algorithm == Scale()
                        iray = -1
                    end

                    # Check for case (2).
                    if iray > 0 && zr + 0.5 * dzr > ztildetfc[ix, jy, kz]

                        # Shift the old ray volume.
                        nray[ix, jy, kz + 1] += 1
                        nrlc = nray[ix, jy, kz + 1]
                        if nrlc > nray_wrk
                            error(
                                "Error in orographic_source: nrlc > nray_wrk!",
                            )
                        end
                        copy_rays!(
                            rays,
                            (iray, ix, jy, kz),
                            (nrlc, ix, jy, kz + 1),
                        )

                        # Clip or extend the old ray volume.
                        if zr - 0.5 * dzr < ztildetfc[ix, jy, kz] || kz2 == 1
                            rays.dzray[nrlc, ix, jy, kz + 1] =
                                zr + 0.5 * dzr - ztildetfc[ix, jy, kz]
                            rays.z[nrlc, ix, jy, kz + 1] =
                                zr + 0.5 * dzr -
                                0.5 * rays.dzray[nrlc, ix, jy, kz + 1]
                            rays.area_zm[nrlc, ix, jy, kz + 1] =
                                rays.dzray[nrlc, ix, jy, kz + 1] *
                                rays.dmray[nrlc, ix, jy, kz + 1]
                        end

                        # Check for case (1).
                    elseif iray < 0
                        nray[ix, jy, kz] += 1
                        iray = nray[ix, jy, kz]
                        if iray > nray_wrk
                            error(
                                "Error in orographic_source: iray > nray_wrk!",
                            )
                        end
                        ir_sfc[i_sfc, ix, jy] = iray
                    end

                    # No ray volume is launched for zero wave-action density.
                else
                    ir_sfc[i_sfc, ix, jy] = -1
                    continue
                end
            end

            # Scale the wave action density.
            if wkb_mode != SteadyState() && launch_algorithm == Scale()
                cgrz =
                    wnrm * (f_cor_nd^2.0 - bvsavg) * wnrh^2.0 / omir /
                    (wnrh^2.0 + wnrm^2.0)^2.0
                wadr *= dt * cgrz / jac[ix, jy, kz] / dz
            end

            # Set physical ray-volume positions.
            rays.x[iray, ix, jy, kz] =
                (x[io + ix] - 0.5 * dx + (ix2 - 0.5) * dx / nrxl)
            rays.y[iray, ix, jy, kz] =
                (y[jo + jy] - 0.5 * dy + (jy2 - 0.5) * dy / nryl)
            rays.z[iray, ix, jy, kz] = (
                ztfc[ix, jy, kz] - 0.5 * jac[ix, jy, kz] * dz +
                (kz2 - 0.5) * jac[ix, jy, kz] * dz / nrzl
            )

            # Set physical ray-volume extent.
            rays.dxray[iray, ix, jy, kz] = dx / nrxl
            rays.dyray[iray, ix, jy, kz] = dy / nryl
            rays.dzray[iray, ix, jy, kz] = jac[ix, jy, kz] * dz / nrzl

            # Compute spectral ray-volume extent.
            if sizex == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd = fac_dk_init * sqrt(wnrk^2 + wnrl^2)
            end
            if sizey == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd = fac_dl_init * sqrt(wnrk^2 + wnrl^2)
            end
            if wnrm == 0.0
                error("Error in orographic_source: wnrm = 0!")
            else
                dm_ini_nd = fac_dm_init * abs[wnrm]
            end

            # Set spectral ray-volume position.
            rays.k[iray, ix, jy, kz] =
                (wnrk - 0.5 * dk_ini_nd + (ik - 0.5) * dk_ini_nd / nrk_init)
            rays.l[iray, ix, jy, kz] =
                (wnrl - 0.5 * dl_ini_nd + (jl - 0.5) * dl_ini_nd / nrl_init)
            rays.m[iray, ix, jy, kz] =
                (wnrm - 0.5 * dm_ini_nd + (km - 0.5) * dm_ini_nd / nrm_init)

            # Set spectral ray-voume extent.
            rays.dkray[iray, ix, jy, kz] = dk_ini_nd / nrk_init
            rays.dlray[iray, ix, jy, kz] = dl_ini_nd / nrl_init
            rays.dmray[iray, ix, jy, kz] = dm_ini_nd / nrm_init

            # Set phase-space volume.
            rays.area_xk[iray, ix, jy, kz] =
                rays.dxray[iray, ix, jy, kz] * rays.dkray[iray, ix, jy, kz]
            rays.area_yl[iray, ix, jy, kz] =
                rays.dyray[iray, ix, jy, kz] * rays.dlray[iray, ix, jy, kz]
            rays.area_zm[iray, ix, jy, kz] =
                rays.dzray[iray, ix, jy, kz] * rays.dmray[iray, ix, jy, kz]

            # Compute spectral volume.
            pspvol = dm_ini_nd
            if sizex > 1
                pspvol *= dk_ini_nd
            end
            if sizey > 1
                pspvol *= dl_ini_nd
            end

            # Set phase-space wave-action density.
            rays.dens[iray, ix, jy, kz] = wadr / pspvol

            # Set intrinsic frequency.
            rays.omega[iray, ix, jy, kz] = omir
        end
    end
end
