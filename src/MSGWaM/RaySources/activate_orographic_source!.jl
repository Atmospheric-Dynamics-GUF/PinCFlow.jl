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

 1. There is a ray volume with nonzero wave-action density that has at least partially passed through the lower boundary. The ray volume is either clipped or extended, such that its lower edge coincides with the surface, and the part below the surface is discarded. Then, it is assigned to the first model layer `k0`, i.e. its indices are changed from `(iray, ix, jy, k0 - 1)` to `(jray, ix, jy, k0)`, where `jray` is the new last ray-volume index at `(ix, jy, k0)`. Finally, a new ray volume is launched.

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
    (; branchr, blocking, long_threshold, nwm) = state.namelists.wkb
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
    @ivy for jy in j0:j1, ix in i0:i1

        # Sum the magnitudes of the spectrum.
        hsum = sum(abs, topography_spectrum[:, ix, jy])

        # Average mean wind, reference density and buoyancy frequency.
        uavg = 0.0
        vavg = 0.0
        rhoavg = 0.0
        bvsavg = 0.0
        dzsum = 0.0
        for kz in k0:k1
            uavg +=
                (u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * jac[ix, jy, kz] * dz
            vavg +=
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * jac[ix, jy, kz] * dz
            rhoavg += rhostrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            bvsavg += bvsstrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            dzsum += jac[ix, jy, kz] * dz
            if ztildetfc[ix, jy, kz] > ztildetfc[ix, jy, k0 - 1] + hsum
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
            zb[ix, jy] = ztildetfc[ix, jy, k0 - 1] + hsum * (1 - 2 * ratio)
        elseif blocking
            ratio = 1
            zb[ix, jy] = ztildetfc[ix, jy, k0 - 1]
        else
            ratio = 1
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
                fc,
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

function activate_orographic_source!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
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
    (; ir_sfc, ix2_sfc, jy2_sfc, kz2_sfc, ik_sfc, jl_sfc, km_sfc, iwm_sfc) =
        state.wkb.surface_indices
    (; nray_wrk, n_sfc, nray, rays, zb, increments) = state.wkb

    if ko != 0
        return
    end

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    # Iterate over surface grid cells.
    @ivy for jy in j0:j1, ix in i0:i1

        # Sum the magnitudes of the spectrum.
        hsum = sum(abs, topography_spectrum[:, ix, jy])

        # Average mean wind, reference density and buoyancy frequency.
        uavg = 0.0
        vavg = 0.0
        rhoavg = 0.0
        bvsavg = 0.0
        dzsum = 0.0
        for kz in k0:k1
            uavg +=
                (u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * jac[ix, jy, kz] * dz
            vavg +=
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * jac[ix, jy, kz] * dz
            rhoavg += rhostrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            bvsavg += bvsstrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            dzsum += jac[ix, jy, kz] * dz
            if ztildetfc[ix, jy, kz] > ztildetfc[ix, jy, k0 - 1] + hsum
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
            zb[ix, jy] = ztildetfc[ix, jy, k0 - 1] + hsum * (1 - 2 * ratio)
        elseif blocking
            ratio = 1.0
            zb[ix, jy] = ztildetfc[ix, jy, k0 - 1]
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
                fc,
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
                if iray < 0
                    if wadr == 0
                        continue
                    else
                        nray[ix, jy, kz] += 1
                        iray = nray[ix, jy, kz]
                        ir_sfc[i_sfc, ix, jy] = iray
                    end
                elseif iray > 0
                    if wadr == 0
                        rays.dens[iray, ix, jy, kz] = 0.0
                        ir_sfc[i_sfc, ix, jy] = -1
                        continue
                    end
                end
            else
                if iray > 0
                    # Shift and clip/extend the old ray volume.
                    if zr + dzr / 2 > ztildetfc[ix, jy, kz]

                        # Shift the old ray volume.
                        nray[ix, jy, kz + 1] += 1
                        nrlc = nray[ix, jy, kz + 1]
                        if nrlc > nray_wrk
                            error(
                                "Error in activate_orographic_source!: nrlc > nray_wrk!",
                            )
                        end
                        copy_rays!(
                            rays,
                            (iray, ix, jy, kz),
                            (nrlc, ix, jy, kz + 1),
                        )
                        for field in fieldnames(WKBIncrements)
                            getfield(increments, field)[nrlc, ix, jy, kz + 1] =
                                getfield(increments, field)[iray, ix, jy, kz]
                            getfield(increments, field)[iray, ix, jy, kz] = 0.0
                        end

                        # Clip/extend the old ray volume.
                        if zr - dzr / 2 < ztildetfc[ix, jy, kz] || kz2 == 1
                            rays.dzray[nrlc, ix, jy, kz + 1] =
                                zr + dzr / 2 - ztildetfc[ix, jy, kz]
                            rays.z[nrlc, ix, jy, kz + 1] =
                                zr + dzr / 2 -
                                rays.dzray[nrlc, ix, jy, kz + 1] / 2
                        end
                    end

                    if wadr == 0
                        rays.dens[iray, ix, jy, kz] = 0.0
                        ir_sfc[i_sfc, ix, jy] = -1
                        continue
                    end
                elseif iray < 0
                    if wadr == 0
                        continue
                    else
                        nray[ix, jy, kz] += 1
                        iray = nray[ix, jy, kz]
                        if iray > nray_wrk
                            error(
                                "Error in activate_orographic_source!: iray > nray_wrk!",
                            )
                        end
                        ir_sfc[i_sfc, ix, jy] = iray
                    end
                end
            end

            # Set physical ray-volume positions.
            rays.x[iray, ix, jy, kz] =
                (x[io + ix] - dx / 2 + (ix2 - 0.5) * dx / nrxl)
            rays.y[iray, ix, jy, kz] =
                (y[jo + jy] - dy / 2 + (jy2 - 0.5) * dy / nryl)
            rays.z[iray, ix, jy, kz] = (
                ztfc[ix, jy, kz] - jac[ix, jy, kz] * dz / 2 +
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
                dm_ini_nd = fac_dm_init * abs(wnrm)
            end

            # Set spectral ray-volume position.
            rays.k[iray, ix, jy, kz] =
                (wnrk - dk_ini_nd / 2 + (ik - 0.5) * dk_ini_nd / nrk_init)
            rays.l[iray, ix, jy, kz] =
                (wnrl - dl_ini_nd / 2 + (jl - 0.5) * dl_ini_nd / nrl_init)
            rays.m[iray, ix, jy, kz] =
                (wnrm - dm_ini_nd / 2 + (km - 0.5) * dm_ini_nd / nrm_init)

            # Set spectral ray-voume extent.
            rays.dkray[iray, ix, jy, kz] = dk_ini_nd / nrk_init
            rays.dlray[iray, ix, jy, kz] = dl_ini_nd / nrl_init
            rays.dmray[iray, ix, jy, kz] = dm_ini_nd / nrm_init

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
        end
    end
end
