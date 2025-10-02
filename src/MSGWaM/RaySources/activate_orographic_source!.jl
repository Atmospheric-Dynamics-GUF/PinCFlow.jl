"""
```julia
activate_orographic_source!(
    state::State,
    omi_ini::AbstractArray{<:AbstractFloat, 4},
    wnk_ini::AbstractArray{<:AbstractFloat, 4},
    wnl_ini::AbstractArray{<:AbstractFloat, 4},
    wnm_ini::AbstractArray{<:AbstractFloat, 4},
    wad_ini::AbstractArray{<:AbstractFloat, 4},
)
```

Compute ray-volume properties in the launch layer (i.e. at `k = k0 - 1`) for the initialization of MSGWaM.

Sets the launch-layer values of arrays for initial ray-volume properties (intrinsic frequencies, wavenumbers and wave-action densities). For this purpose, the horizontal components of the resolved wind, the background density and the squared buoyancy frequency are vertically averaged between the surface and an approximation for the summits of the unresolved orography. The vertical averages are then used to compute a non-dimensionalized mountain wave amplitude, from which an approximate reduction of the generated wave amplitude due to blocking is inferred (see below). Afterwards, the ray-volume properties are obtained by calling `compute_orographic_mode` with the correspondingly scaled mode of the orographic spectrum and the vertical averages as arguments.

```julia
activate_orographic_source!(state::State)
```

Launch ray volumes that represent unresolved orographic gravity waves.

In each column of MPI processes at the lower boundary, this method first computes vertical averages of the horizontal components of the resolved wind, the background density and the squared buoyancy frequency between ``h_\\mathrm{b}`` (the surface) and ``h_\\mathrm{b} + \\Delta h`` (an approximation for the summits of the unresolved orography, with ``\\Delta h = \\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right|``). The vertical averages are then used to compute a non-dimensionalized mountain wave amplitude, from which an approximate reduction of the generated wave amplitude due to blocking, as well as the depth of the blocked layer, is inferred. A loop over the spectral modes of the unresolved orography follows, in which the properties of each mode are computed, using `compute_orographic_mode` with the scaled mode of the orographic spectrum and vertical averages as arguments, and corresponding ray volumes are launched at `k = k0 - 1`.

The parameterization of blocking is built around the non-dimensionalized mountain wave amplitude, or Long number,

```math
\\mathrm{Lo} = \\frac{N_h \\Delta h}{\\left|\\boldsymbol{u}_h\\right|},
```

where ``N_h`` is the square root of the vertically averaged squared buoyancy frequency and ``\\boldsymbol{u}_h`` is the vertically averaged resolved horizontal wind. This number is used to estimate the depth of the blocked layer as

```math
\\Delta z_\\mathrm{B} = 2 \\Delta h \\max \\left(0, \\frac{\\mathrm{Lo} - C}{\\mathrm{Lo}}\\right),
```

where ``C`` is a critical value represented by the model parameter `state.namelists.wkb.long_threshold`. The corresponding scaling of the orographic spectrum is given by

```math
r \\left(\\mathrm{Lo}\\right) = \\frac{2 \\Delta h - \\Delta z_\\mathrm{B}}{2 \\Delta h} = \\min \\left(1, \\frac{C}{\\mathrm{Lo}}\\right),
```

so that ``\\Delta z_\\mathrm{B} = 2 \\Delta h \\left(1 - r\\right)``. In addition to the reduction of the mountain-wave amplitude, the present blocked-layer scheme adds a blocked-flow drag to the mean-flow impact. This is implemented in [`PinCFlow.MSGWaM.MeanFlowEffect.apply_blocked_layer_scheme!`](@ref).

The launch algorithm distinguishes between the following situations (regarding previously launched ray volumes).

 1. There is no ray volume with nonzero wave-action density. A new ray volume is launched.

 1. There is a ray volume with nonzero wave-action density that has at least partially passed through the lower boundary. The ray volume is either clipped or extended, such that its lower edge coincides with the surface, and the part below the surface is discarded. Then, it is assigned to the first model layer `k0`, i.e. its indices are changed from `(r, i, j, k0 - 1)` to `(rray, i, j, k0)`, where `rray` is the new last ray-volume index at `(i, j, k0)`. Finally, a new ray volume is launched.

 1. There is a ray volume with nonzero wave-action density, which has not yet crossed the lower boundary. It is replaced with a new one.

# Arguments

  - `state`: Model state.

  - `omi_ini`: Array for intrinsic frequencies.

  - `wnk_ini`: Array for zonal wavenumbers.

  - `wnl_ini`: Array for meridional wavenumbers.

  - `wnm_ini`: Array for vertical wavenumbers.

  - `wad_ini`: Array for wave-action densities.

# See also

  - [`PinCFlow.MSGWaM.RaySources.compute_orographic_mode`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function activate_orographic_source! end

function activate_orographic_source!(
    state::State,
    omi_ini::AbstractArray{<:AbstractFloat, 4},
    wnk_ini::AbstractArray{<:AbstractFloat, 4},
    wnl_ini::AbstractArray{<:AbstractFloat, 4},
    wnm_ini::AbstractArray{<:AbstractFloat, 4},
    wad_ini::AbstractArray{<:AbstractFloat, 4},
)
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch, blocking, long_threshold, wave_modes) = state.namelists.wkb
    (; tref) = state.constants
    (; ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, ztildetfc, k_spectrum, l_spectrum, topography_spectrum) =
        state.grid
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; u, v) = state.variables.predictands
    (; zb) = state.wkb

    if ko != 0
        return
    end

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    # Iterate over surface grid cells.
    @ivy for j in j0:j1, i in i0:i1

        # Sum the magnitudes of the spectrum.
        hsum = sum(abs, topography_spectrum[:, i, j])

        # Average mean wind, reference density and buoyancy frequency.
        uavg = 0.0
        vavg = 0.0
        rhoavg = 0.0
        bvsavg = 0.0
        dzsum = 0.0
        for k in k0:k1
            uavg += (u[i, j, k] + u[i - 1, j, k]) / 2 * jac[i, j, k] * dz
            vavg += (v[i, j, k] + v[i, j - 1, k]) / 2 * jac[i, j, k] * dz
            rhoavg += rhostrattfc[i, j, k] * jac[i, j, k] * dz
            bvsavg += bvsstrattfc[i, j, k] * jac[i, j, k] * dz
            dzsum += jac[i, j, k] * dz
            if ztildetfc[i, j, k] > ztildetfc[i, j, k0 - 1] + hsum
                break
            end
        end
        uavg /= dzsum
        vavg /= dzsum
        rhoavg /= dzsum
        bvsavg /= dzsum

        # Determine the blocked layer.
        if blocking && hsum > 0
            long = sqrt(bvsavg) / sqrt(uavg^2 + vavg^2) * hsum
            ratio = min(1, long_threshold / long)
            zb[i, j] = ztildetfc[i, j, k0 - 1] + hsum * (1 - 2 * ratio)
        elseif blocking
            ratio = 1
            zb[i, j] = ztildetfc[i, j, k0 - 1]
        else
            ratio = 1
        end

        # Set launch level.
        k = k0 - 1

        # Iterate over wave modes.
        for alpha in 1:wave_modes

            # Compute intrinsic frequency, wavenumbers and wave-action density.
            (omi, wnk, wnl, wnm, wad) = compute_orographic_mode(
                ratio * topography_spectrum[alpha, i, j],
                k_spectrum[alpha, i, j],
                l_spectrum[alpha, i, j],
                uavg,
                vavg,
                rhoavg,
                bvsavg,
                fc,
                branch,
            )

            # Save the results.
            omi_ini[alpha, i, j, k] = omi
            wnk_ini[alpha, i, j, k] = wnk
            wnl_ini[alpha, i, j, k] = wnl
            wnm_ini[alpha, i, j, k] = wnm
            wad_ini[alpha, i, j, k] = wad
        end
    end
    return
end

function activate_orographic_source!(state::State)
    (; ndx, ndy) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (;
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        fdk,
        fdl,
        fdm,
        branch,
        blocking,
        long_threshold,
        wkb_mode,
    ) = state.namelists.wkb
    (; tref) = state.constants
    (; io, jo, ko, i0, i1, j0, j1, k0, k1) = state.domain
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
    (; rs, ixs, jys, kzs, iks, jls, kms, alphas) = state.wkb.surface_indices
    (; nray_wrk, n_sfc, nray, rays, zb, increments) = state.wkb

    if ko != 0
        return
    end

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    # Iterate over surface grid cells.
    @ivy for j in j0:j1, i in i0:i1

        # Sum the magnitudes of the spectrum.
        hsum = sum(abs, topography_spectrum[:, i, j])

        # Average mean wind, reference density and buoyancy frequency.
        uavg = 0.0
        vavg = 0.0
        rhoavg = 0.0
        bvsavg = 0.0
        dzsum = 0.0
        for k in k0:k1
            uavg += (u[i, j, k] + u[i - 1, j, k]) / 2 * jac[i, j, k] * dz
            vavg += (v[i, j, k] + v[i, j - 1, k]) / 2 * jac[i, j, k] * dz
            rhoavg += rhostrattfc[i, j, k] * jac[i, j, k] * dz
            bvsavg += bvsstrattfc[i, j, k] * jac[i, j, k] * dz
            dzsum += jac[i, j, k] * dz
            if ztildetfc[i, j, k] > ztildetfc[i, j, k0 - 1] + hsum
                break
            end
        end
        uavg /= dzsum
        vavg /= dzsum
        rhoavg /= dzsum
        bvsavg /= dzsum

        # Determine the blocked layer.
        if blocking && hsum > 0
            long = sqrt(bvsavg) / sqrt(uavg^2 + vavg^2) * hsum
            ratio = min(1, long_threshold / long)
            zb[i, j] = ztildetfc[i, j, k0 - 1] + hsum * (1 - 2 * ratio)
        elseif blocking
            ratio = 1.0
            zb[i, j] = ztildetfc[i, j, k0 - 1]
        else
            ratio = 1.0
        end

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
            (omir, wnrk, wnrl, wnrm, wadr) = compute_orographic_mode(
                ratio * topography_spectrum[alpha, i, j],
                k_spectrum[alpha, i, j],
                l_spectrum[alpha, i, j],
                uavg,
                vavg,
                rhoavg,
                bvsavg,
                fc,
                branch,
            )

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
            if wkb_mode == SteadyState()
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
                    if zr + dzr / 2 > ztildetfc[i, j, k]

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
                        for field in fieldnames(WKBIncrements)
                            getfield(increments, field)[
                                local_count,
                                i,
                                j,
                                k + 1,
                            ] = getfield(increments, field)[r, i, j, k]
                            getfield(increments, field)[r, i, j, k] = 0.0
                        end

                        # Clip/extend the old ray volume.
                        if zr - dzr / 2 < ztildetfc[i, j, k] || kz == 1
                            rays.dzray[local_count, i, j, k + 1] =
                                zr + dzr / 2 - ztildetfc[i, j, k]
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
            rays.x[r, i, j, k] = (x[io + i] - dx / 2 + (ix - 0.5) * dx / nrx)
            rays.y[r, i, j, k] = (y[jo + j] - dy / 2 + (jy - 0.5) * dy / nry)
            rays.z[r, i, j, k] = (
                ztfc[i, j, k] - jac[i, j, k] * dz / 2 +
                (kz - 0.5) * jac[i, j, k] * dz / nrz
            )

            # Set physical ray-volume extent.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

            # Compute spectral ray-volume extent.
            if ndx == 1
                dk_ini_nd = 0.0
            else
                dk_ini_nd = fdk * sqrt(wnrk^2 + wnrl^2)
            end
            if ndy == 1
                dl_ini_nd = 0.0
            else
                dl_ini_nd = fdl * sqrt(wnrk^2 + wnrl^2)
            end
            if wnrm == 0.0
                error("Error in orographic_source: wnrm = 0!")
            else
                dm_ini_nd = fdm * abs(wnrm)
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
            if ndx > 1
                pspvol *= dk_ini_nd
            end
            if ndy > 1
                pspvol *= dl_ini_nd
            end

            # Set phase-space wave-action density.
            rays.dens[r, i, j, k] = wadr / pspvol
        end
    end
end
