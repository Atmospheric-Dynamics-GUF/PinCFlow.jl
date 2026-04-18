"""
```julia
activate_orographic_source!(state::State)
```

Launch ray volumes that represent unresolved orographic gravity waves.

This method first calls `compute_orographic_modes!` and then launches corresponding ray volumes. The launch algorithm distinguishes between the following situations (regarding previously launched ray volumes).

 1. There is no ray volume with nonzero wave-action density. A new ray volume is launched.

 1. There is a ray volume with nonzero wave-action density that has at least partially passed through the lower boundary. The ray volume is either clipped or extended, such that its lower edge coincides with the surface, and the part below the surface is discarded. Then, it is assigned to the first model layer `k0`, i.e. its indices are changed from `(r, i, j, k0 - 1)` to `(rray, i, j, k0)`, where `rray` is the new last ray-volume index at `(i, j, k0)`. Finally, a new ray volume is launched.

 1. There is a ray volume with nonzero wave-action density, which has not yet crossed the lower boundary. It is replaced with a new one.

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.MSGWaM.RaySources.compute_orographic_modes!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_increments!`](@ref)

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been validated yet.
"""
function activate_orographic_source! end

function activate_orographic_source!(state::State)
    (; x_size, y_size) = state.namelists.domain
    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        wkb_mode,
    ) = state.namelists.wkb
    (; ko, i0, i1, j0, j1, k0) = state.domain
    (; dx, dy, dz, x, y, zc, jac, zctilde) = state.grid
    (; rs, ixs, jys, kzs, iks, jls, kms, alphas) = state.wkb.surface_indices
    (; nray_wrk, n_sfc, nray, rays, increments, spectrum) = state.wkb

    if ko != 0
        return
    end

    compute_orographic_modes!(state)

    # Loop over surface grid cells.
    @ivy for j in j0:j1, i in i0:i1

        # Set launch level.
        k = k0 - 1

        # Loop over surface ray volumes.
        for s in 1:n_sfc
            r = rs[s, i, j]

            # Set surface indices.
            ix = ixs[s]
            jy = jys[s]
            kz = kzs[s]
            ik = iks[s]
            jl = jls[s]
            km = kms[s]
            alpha = alphas[s]

            # Compute intrinsic frequency, wavenumbers and wave-action density.
            wnrk = spectrum.k[alpha, i, j, k]
            wnrl = spectrum.l[alpha, i, j, k]
            wnrm = spectrum.m[alpha, i, j, k]
            omir = spectrum.omega[alpha, i, j, k]
            wadr = spectrum.a[alpha, i, j, k]

            # Get vertical position and extent of old ray volume.
            if r > 0
                zr = rays.z[r, i, j, k]
                dzr = rays.dzray[r, i, j, k]
            else
                zr = 0.0
                dzr = 0.0
            end

            # Check if a new ray volume is to be launched and clip/extend the
            # old one. Three cases are distinguished.
            # (1) There is no ray volume with nonzero wave-action density. A new
            #     ray volume is launched.
            # (2) There is a ray volume with nonzero wave-action density, which
            #     has partially passed the lower boundary. It is clipped and the
            #     part below the lower boundary discarded before a new ray
            #     volume is launched.
            # (3) There is a ray volume with nonzero wave-action density, which
            #     has not yet crossed the lower boundary. It is replaced with a
            #     new one.
            if wkb_mode == :SteadyState
                if r < 0
                    if wadr == 0
                        continue
                    else
                        nray[i, j, k] += 1
                        r = nray[i, j, k]
                        rs[s, i, j] = r
                    end
                elseif r > 0
                    if wadr == 0
                        rays.dens[r, i, j, k] = 0.0
                        rs[s, i, j] = -1
                        continue
                    end
                end
            else
                if r > 0
                    # Shift and clip/extend the old ray volume.
                    if zr + dzr / 2 > zctilde[i, j, k]

                        # Shift the old ray volume.
                        nray[i, j, k + 1] += 1
                        local_count = nray[i, j, k + 1]
                        if local_count > nray_wrk
                            error(
                                "Error in activate_orographic_source!: local_count > nray_wrk!",
                            )
                        end
                        copy_rays!(
                            rays,
                            r => local_count,
                            i => i,
                            j => j,
                            k => k + 1,
                        )
                        copy_increments!(
                            increments,
                            r => local_count,
                            i => i,
                            j => j,
                            k => k + 1,
                        )
                        for field in fieldnames(WKBIncrements)
                            getfield(increments, field)[r, i, j, k] = 0.0
                        end

                        # Clip/extend the old ray volume.
                        if zr - dzr / 2 < zctilde[i, j, k] || kz == 1
                            rays.dzray[local_count, i, j, k + 1] =
                                zr + dzr / 2 - zctilde[i, j, k]
                            rays.z[local_count, i, j, k + 1] =
                                zr + dzr / 2 -
                                rays.dzray[local_count, i, j, k + 1] / 2
                        end
                    end

                    if wadr == 0
                        rays.dens[r, i, j, k] = 0.0
                        rs[s, i, j] = -1
                        continue
                    end
                elseif r < 0
                    if wadr == 0
                        continue
                    else
                        nray[i, j, k] += 1
                        r = nray[i, j, k]
                        if r > nray_wrk
                            error(
                                "Error in activate_orographic_source!: r > nray_wrk!",
                            )
                        end
                        rs[s, i, j] = r
                    end
                end
            end

            # Set physical ray-volume positions.
            rays.x[r, i, j, k] = (x[i] - dx / 2 + (ix - 0.5) * dx / nrx)
            rays.y[r, i, j, k] = (y[j] - dy / 2 + (jy - 0.5) * dy / nry)
            rays.z[r, i, j, k] = (
                zc[i, j, k] - jac[i, j, k] * dz / 2 +
                (kz - 0.5) * jac[i, j, k] * dz / nrz
            )

            # Set physical ray-volume extent.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

            # Compute spectral ray-volume extent.
            if x_size == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd = dkr_factor * sqrt(wnrk^2 + wnrl^2)
            end
            if y_size == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd = dlr_factor * sqrt(wnrk^2 + wnrl^2)
            end
            if wnrm == 0.0
                error("Error in orographic_source: wnrm = 0!")
            else
                dm_ini_nd = dmr_factor * abs(wnrm)
            end

            # Set spectral ray-volume position.
            rays.k[r, i, j, k] =
                (wnrk - dk_ini_nd / 2 + (ik - 0.5) * dk_ini_nd / nrk)
            rays.l[r, i, j, k] =
                (wnrl - dl_ini_nd / 2 + (jl - 0.5) * dl_ini_nd / nrl)
            rays.m[r, i, j, k] =
                (wnrm - dm_ini_nd / 2 + (km - 0.5) * dm_ini_nd / nrm)

            # Set spectral ray-volume extent.
            rays.dkray[r, i, j, k] = dk_ini_nd / nrk
            rays.dlray[r, i, j, k] = dl_ini_nd / nrl
            rays.dmray[r, i, j, k] = dm_ini_nd / nrm

            # Compute spectral volume.
            pspvol = dm_ini_nd
            if x_size > 1
                pspvol *= dk_ini_nd
            end
            if y_size > 1
                pspvol *= dl_ini_nd
            end

            # Set phase-space wave-action density.
            rays.dens[r, i, j, k] = wadr / pspvol
        end
    end

    return
end
