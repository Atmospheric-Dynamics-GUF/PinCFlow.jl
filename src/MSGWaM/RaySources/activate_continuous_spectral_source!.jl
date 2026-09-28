"""
    activate_continuous_spectral_source!(state::State)

Activate the continuously forced lower-boundary spectral source.

Source ray volumes are stored in the artificial lower layer `k0 - 1`.
After a complete WKB propagation step, every source ray whose upper edge
has crossed the lower physical boundary is copied into the first physical
layer `k0`. Only the portion lying above the boundary is retained.

After the transfer test, every source slot is restored to its prescribed
initial position, spectral properties, extents, and Gaussian wave-action
density. Thus the artificial layer acts as a stationary spectral reservoir.

The vertical source discretization supports arbitrary `nrz >= 1`.
"""
function activate_continuous_spectral_source! end

function activate_continuous_spectral_source!(state::State)
    (; x_size, y_size) = state.namelists.domain

    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        initial_wave_field,
        source_mode,
    ) = state.namelists.wkb

    (; m_sigma_cutoff) = state.namelists.triad
    (; lref, tref, rhoref, uref) = state.constants

    (; ko, i0, i1, j0, j1, k0) = state.domain
    (; dx, dy, dz, x, y, zc, zctilde, jac) = state.grid

    (; nray_wrk, n_sfc, nray, rays, surface_indices, increments) = state.wkb
    (; rs, ixs, jys, kzs, iks, jls, kms, alphas) = surface_indices

    # This routine applies only to the prescribed continuous source.
    if !(source_mode isa ContinuousSpectralSource)
        return
    end

    # The source exists only on the lowest vertical MPI subdomain.
    if ko != 0
        return
    end

    if m_sigma_cutoff <= 0.0
        error(
            "Error in activate_continuous_spectral_source!: ",
            "m_sigma_cutoff must be positive.",
        )
    end

    # ------------------------------------------------------------------
    # Gaussian weights in m.
    #
    # This must be identical to the normalization used during source
    # initialization:
    #
    #     sum(m_weights) / nrm = 1.
    # ------------------------------------------------------------------

    m_weights = zeros(nrm)

    for km in 1:nrm
        xi =
            -m_sigma_cutoff +
            (km - 0.5) * 2.0 * m_sigma_cutoff / nrm

        m_weights[km] = exp(-0.5 * xi^2)
    end

    m_weights .*= nrm / sum(m_weights)

    # Artificial source layer and its upper boundary.
    k = k0 - 1

    # ------------------------------------------------------------------
    # Iterate over horizontal source cells.
    # ------------------------------------------------------------------

    @ivy for j in j0:j1, i in i0:i1

        z_boundary = zctilde[i, j, k]

        # --------------------------------------------------------------
        # Process every persistent source slot.
        # --------------------------------------------------------------

        for s in 1:n_sfc

            r = rs[s, i, j]

            ix = ixs[s]
            jy = jys[s]
            kz = kzs[s]
            ik = iks[s]
            jl = jls[s]
            km = kms[s]
            alpha = alphas[s]

            # ----------------------------------------------------------
            # Reconstruct the prescribed carrier spectrum at the lower
            # physical boundary.
            # ----------------------------------------------------------

            (kdim, ldim, mdim, _, adim) = initial_wave_field(
                alpha,
                x[i] * lref,
                y[j] * lref,
                z_boundary * lref,
            )

            wnk0 = kdim * lref
            wnl0 = ldim * lref
            wnm0 = mdim * lref
            wad0 = adim / rhoref / uref^2 / tref

            # ----------------------------------------------------------
            # If the prescribed source is zero, deactivate this slot.
            # ----------------------------------------------------------

            if wad0 == 0.0
                if r > 0
                    rays.dens[r, i, j, k] = 0.0

                    for field in fieldnames(WKBIncrements)
                        getfield(increments, field)[r, i, j, k] = 0.0
                    end
                end

                rs[s, i, j] = -1
                continue
            end

            # ----------------------------------------------------------
            # If a source ray already exists, check whether any part of
            # it has crossed the lower physical boundary.
            # ----------------------------------------------------------

            if r > 0
                zr = rays.z[r, i, j, k]
                dzr = rays.dzray[r, i, j, k]

                z_lower = zr - dzr / 2
                z_upper = zr + dzr / 2

                if z_upper > z_boundary

                    # --------------------------------------------------
                    # Create a new physical ray in the first model layer.
                    # --------------------------------------------------

                    nray[i, j, k0] += 1
                    rnew = nray[i, j, k0]

                    if rnew > nray_wrk
                        error(
                            "Error in activate_continuous_spectral_source!: ",
                            "Number of ray volumes in first physical layer ",
                            "exceeds nray_wrk = ",
                            nray_wrk,
                            " at ",
                            (i, j, k0),
                        )
                    end

                    copy_rays!(
                        rays,
                        r => rnew,
                        i => i,
                        j => j,
                        k => k0,
                    )

                    # --------------------------------------------------
                    # Retain only the part of the ray that has entered
                    # the physical domain.
                    #
                    # If the lower edge is still below the boundary,
                    # clip it at the boundary.
                    #
                    # For the lowest source subvolume kz == 1, also
                    # extend back to the boundary if the whole source
                    # layer has crossed during one time step. This is
                    # the same gap-prevention logic used by the existing
                    # orographic source.
                    # --------------------------------------------------

                    if z_lower < z_boundary || kz == 1
                        dz_phys = z_upper - z_boundary

                        if dz_phys <= 0.0
                            error(
                                "Error in activate_continuous_spectral_source!: ",
                                "Non-positive transferred vertical extent.",
                            )
                        end

                        rays.dzray[rnew, i, j, k0] = dz_phys
                        rays.z[rnew, i, j, k0] =
                            z_boundary + dz_phys / 2
                    end

                    # The physical ray starts a new WKB step with zero
                    # Runge-Kutta history. The next rkstage == 1 would
                    # reset this anyway, but doing it explicitly avoids
                    # carrying source-layer increments into k0.
                    for field in fieldnames(WKBIncrements)
                        getfield(increments, field)[rnew, i, j, k0] = 0.0
                    end
                end
            else
                # ------------------------------------------------------
                # Recreate an absent source slot.
                # ------------------------------------------------------

                nray[i, j, k] += 1
                r = nray[i, j, k]

                if r > nray_wrk
                    error(
                        "Error in activate_continuous_spectral_source!: ",
                        "Number of source ray volumes exceeds nray_wrk = ",
                        nray_wrk,
                        " at ",
                        (i, j, k),
                    )
                end

                rs[s, i, j] = r
            end

            # ----------------------------------------------------------
            # Restore the source slot.
            #
            # This happens whether or not the previous source ray crossed.
            # Thus every WKB time step begins with the same prescribed
            # source reservoir.
            # ----------------------------------------------------------

            rays.x[r, i, j, k] =
                x[i] -
                0.5 * dx +
                (ix - 0.5) * dx / nrx

            rays.y[r, i, j, k] =
                y[j] -
                0.5 * dy +
                (jy - 0.5) * dy / nry

            rays.z[r, i, j, k] =
                zc[i, j, k] -
                0.5 * jac[i, j, k] * dz +
                (kz - 0.5) * jac[i, j, k] * dz / nrz

            # Physical extents.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] =
                jac[i, j, k] * dz / nrz

            # ----------------------------------------------------------
            # Reconstruct spectral widths.
            # ----------------------------------------------------------

            wnh0 = sqrt(wnk0^2 + wnl0^2)

            if wnh0 <= 0.0
                error(
                    "Error in activate_continuous_spectral_source!: ",
                    "Horizontal source wavenumber must be nonzero ",
                    "for mode ",
                    alpha,
                    ".",
                )
            end

            if wnm0 == 0.0
                error(
                    "Error in activate_continuous_spectral_source!: ",
                    "Source vertical wavenumber is zero for mode ",
                    alpha,
                    ".",
                )
            end

            if x_size == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd =
                    dkr_factor[alpha] * wnh0

                if dk_ini_nd <= 0.0
                    error(
                        "Error in activate_continuous_spectral_source!: ",
                        "dk_ini_nd <= 0 for source mode ",
                        alpha,
                        ".",
                    )
                end
            end

            if y_size == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd =
                    dlr_factor[alpha] * wnh0

                if dl_ini_nd <= 0.0
                    error(
                        "Error in activate_continuous_spectral_source!: ",
                        "dl_ini_nd <= 0 for source mode ",
                        alpha,
                        ".",
                    )
                end
            end

            dm_ini_nd =
                dmr_factor[alpha] * abs(wnm0)

            if dm_ini_nd <= 0.0
                error(
                    "Error in activate_continuous_spectral_source!: ",
                    "dm_ini_nd <= 0 for source mode ",
                    alpha,
                    ".",
                )
            end

            # ----------------------------------------------------------
            # Restore prescribed spectral positions and extents.
            # ----------------------------------------------------------

            rays.k[r, i, j, k] =
                wnk0 -
                0.5 * dk_ini_nd +
                (ik - 0.5) * dk_ini_nd / nrk

            rays.l[r, i, j, k] =
                wnl0 -
                0.5 * dl_ini_nd +
                (jl - 0.5) * dl_ini_nd / nrl

            rays.m[r, i, j, k] =
                wnm0 -
                0.5 * dm_ini_nd +
                (km - 0.5) * dm_ini_nd / nrm

            rays.dkray[r, i, j, k] =
                dk_ini_nd / nrk

            rays.dlray[r, i, j, k] =
                dl_ini_nd / nrl

            rays.dmray[r, i, j, k] =
                dm_ini_nd / nrm

            # Complete spectral volume.
            pspvol = dm_ini_nd

            if x_size > 1
                pspvol *= dk_ini_nd
            end

            if y_size > 1
                pspvol *= dl_ini_nd
            end

            # Restore the prescribed Gaussian wave-action density.
            rays.dens[r, i, j, k] =
                wad0 / pspvol * m_weights[km]

            # The restored reservoir starts the next WKB step without
            # accumulated Runge-Kutta increments.
            for field in fieldnames(WKBIncrements)
                getfield(increments, field)[r, i, j, k] = 0.0
            end
        end
    end

    return
end