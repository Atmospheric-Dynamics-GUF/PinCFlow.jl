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

Initialize orographic wave source for ray initialization phase.

Computes initial wave mode characteristics for all wavenumber modes at the
surface level, used during ray volume initialization.

# Arguments

  - `state::State`: Complete simulation state
  - `omi_ini`: Output array for intrinsic frequencies [nwm, nx, ny, nz]
  - `wnk_ini`: Output array for x-direction wavenumbers [nwm, nx, ny, nz]
  - `wnl_ini`: Output array for y-direction wavenumbers [nwm, nx, ny, nz]
  - `wnm_ini`: Output array for vertical wavenumbers [nwm, nx, ny, nz]
  - `wad_ini`: Output array for wave action densities [nwm, nx, ny, nz]

# Process

 1. **Surface Only**: Only processes surface grid cells (`ko == 0`)
 2. **Background Averaging**: Computes column-averaged wind, density, and stratification
 3. **Topographic Spectrum**: Uses pre-computed topographic wavenumber spectrum
 4. **Mode Calculation**: Calls `compute_orographic_mode` for each spectral mode
 5. **Storage**: Saves results in provided arrays for later ray initialization

# Applications

Used during simulation initialization to pre-compute source characteristics.
"""
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
                (u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * jac[ix, jy, kz] * dz
            vavg +=
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * jac[ix, jy, kz] * dz
            rhoavg += rhostrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            bvsavg += bvsstrattfc[ix, jy, kz] * jac[ix, jy, kz] * dz
            dzsum += jac[ix, jy, kz] * dz
            @views if ztildetfc[ix, jy, kz] >
                      ztildetfc[ix, jy, k0 - 1] +
                      sum(abs.(topography_spectrum[:, ix, jy]))
                break
            end
        end
        uavg /= dzsum
        vavg /= dzsum
        rhoavg /= dzsum
        bvsavg /= dzsum

        # Determine the blocked layer.
        @views if blocking && sum(abs.(topography_spectrum[:, ix, jy])) > 0
            long =
                sqrt(bvsavg) / sqrt(uavg^2 + vavg^2) *
                sum(abs.(topography_spectrum[:, ix, jy]))
            ratio = min(1, long_threshold / long)
            zb[ix, jy] =
                ztildetfc[ix, jy, k0 - 1] +
                sum(abs.(topography_spectrum[:, ix, jy])) * (1 - 2 * ratio)
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

"""
```julia
activate_orographic_source!(state::State)
```

Activate orographic wave source during time integration.

Updates surface ray volumes to maintain continuous wave generation from
topography, handling ray replacement, clipping, and new ray creation
as rays propagate away from the surface.

# Arguments

  - `state::State`: Complete simulation state

# Algorithm

 1. **Surface Processing**: Only operates on surface processes (`ko == 0`)

 2. **Background Conditions**: Re-computes column-averaged atmospheric state
 3. **Blocking Calculation**: Applies flow blocking parameterization if enabled
 4. **Ray Management**: For each surface ray volume:

      + **Existing rays**: Update or remove based on wave conditions
      + **Propagated rays**: Clip/extend rays crossing vertical boundaries
      + **New rays**: Create fresh ray volumes with current source strength
 5. **Mode Characteristics**: Re-compute wave properties for current conditions
 6. **Ray Properties**: Set physical/spectral positions, extents, and wave action

# Flow Blocking

When enabled, reduces effective topographic height based on:

  - **Linearity parameter**: `L = N·h / U` (ratio of buoyancy to advection timescales)
  - **Blocking ratio**: `r = min(1, L_threshold / L)`
  - **Effective height**: `h_eff = h · r`
  - **Blocked layer depth**: `z_b = h · (1 - 2r)`

# Ray Volume Management

Three cases for existing ray volumes:

 1. **No existing ray** (`iray < 0`): Create new ray if source is active
 2. **Existing ray, zero source**: Remove ray and mark inactive
 3. **Existing ray, active source**: Update ray properties

# Boundary Handling

  - **Steady State**: Simple replacement of ray properties

  - **Time Dependent**:

      + Shift rays crossing vertical boundaries to next level
      + Clip ray extents at boundary interfaces
      + Preserve wave action through coordinate transformations

# Wave Action Conservation

Maintains consistency between:

  - Source strength from topographic spectrum
  - Ray volume phase space density
  - Spectral volume normalization factors

# Error Checking

  - Validates ray counts don't exceed working limits
  - Ensures non-zero vertical wavenumber for active sources
  - Checks ray positioning within valid domains
"""
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
    for jy in j0:j1, ix in i0:i1

        # Average mean wind, reference density and buoyancy frequency. This
        # should be done without a vertical loop.
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
            @views if ztildetfc[ix, jy, kz] >
                      ztildetfc[ix, jy, k0 - 1] +
                      sum(abs.(topography_spectrum[:, ix, jy]))
                break
            end
        end
        uavg /= dzsum
        vavg /= dzsum
        rhoavg /= dzsum
        bvsavg /= dzsum

        # Determine the blocked layer.
        @views if blocking && sum(abs.(topography_spectrum[:, ix, jy])) > 0
            long =
                sqrt(bvsavg) / sqrt(uavg^2 + vavg^2) *
                sum(abs.(topography_spectrum[:, ix, jy]))
            ratio = min(1, long_threshold / long)
            zb[ix, jy] =
                ztildetfc[ix, jy, k0 - 1] +
                sum(abs.(topography_spectrum[:, ix, jy])) * (1 - 2 * ratio)
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
                        for field in fieldnames(Increments)
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
