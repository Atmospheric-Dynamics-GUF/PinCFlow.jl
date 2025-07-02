"""
    compute_gw_integrals!(state::State, wkb_mode::MultiColumn)

Compute gravity wave momentum flux integrals for multi-column WKB mode.

Calculates momentum flux and Eliassen-Palm flux components by integrating
over all ray volumes, accounting for their spatial and spectral distribution
across grid cells.

# Arguments

  - `state::State`: Complete simulation state
  - `wkb_mode::MultiColumn`: Multi-column WKB mode allowing horizontal propagation

# Computed Integrals

## Momentum Flux Components

  - `uu`: Zonal momentum flux in x-direction `⟨u'u'⟩`
  - `uv`: Cross momentum flux `⟨u'v'⟩`
  - `uw`: Vertical zonal momentum flux `⟨u'w'⟩`
  - `vv`: Meridional momentum flux in y-direction `⟨v'v'⟩`
  - `vw`: Vertical meridional momentum flux `⟨v'w'⟩`

## Eliassen-Palm Components

  - `etx`: Zonal EP flux component
  - `ety`: Meridional EP flux component
  - `utheta`: Zonal heat flux `⟨u'θ'⟩`
  - `vtheta`: Meridional heat flux `⟨v'θ'⟩`

## Energy

  - `e`: Wave energy density

# Integration Process

For each ray volume:

 1. **Overlap Calculation**: Determine which grid cells the ray intersects
 2. **Phase Space Factors**: Compute spectral and spatial overlap fractions
 3. **Group Velocities**: Calculate from local dispersion relation
 4. **Flux Components**: Compute momentum flux using group velocity relations
 5. **Grid Integration**: Distribute ray contribution across overlapping cells

# Dispersion Relations

Uses full gravity wave dispersion relation including rotation:

  - `ω² = N²k_h²/(k_h² + m²) + f²m²/(k_h² + m²)`
  - Group velocities: `cg_i = ∂ω/∂k_i`

# Rotational Effects

When Coriolis parameter `f ≠ 0`:

  - Modified momentum flux expressions account for wave polarization
  - Eliassen-Palm fluxes include geostrophic adjustment effects
  - Heat flux components computed from EP flux relations

# Spatial Distribution

Ray contributions are distributed across grid cells based on:

  - **Physical overlap**: Ray volume intersection with grid cells
  - **Spectral weighting**: Wavenumber extent and resolution
  - **Phase space volume**: Conservation of wave action density    # Set Coriolis parameter.
"""
function compute_gw_integrals!(state::State, wkb_mode::MultiColumn)
    (; domain, grid) = state
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = domain
    (; dx, dy, dz, x, y, ztildetfc, jac) = grid
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(GWIntegrals)
        getfield(integrals, field) .= 0.0
    end

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            (ixmin, ixmax, jymin, jymax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + ix] + dx / 2) -
                        max(xr - dxr / 2, x[io + ix] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jy] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jy] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        if sizex > 1
                            if fc != 0
                                integrals.uu[ix, jy, kz] +=
                                    wadr * (
                                        kr * cgirx -
                                        (kr * cgirx + lr * cgiry) /
                                        (1 - (omir / fc)^2)
                                    )
                            else
                                integrals.uu[ix, jy, kz] += wadr * kr * cgirx
                            end
                        end

                        if sizex > 1 || sizey > 1
                            integrals.uv[ix, jy, kz] += wadr * cgirx * lr
                        end

                        integrals.uw[ix, jy, kz] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        if sizey > 1
                            if fc != 0
                                integrals.vv[ix, jy, kz] +=
                                    wadr * (
                                        lr * cgiry -
                                        (kr * cgirx + lr * cgiry) /
                                        (1 - (omir / fc)^2)
                                    )
                            else
                                integrals.vv[ix, jy, kz] += wadr * lr * cgiry
                            end
                        end

                        integrals.vw[ix, jy, kz] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.etx[ix, jy, kz] +=
                                wadr * fc^2 * n2r * kr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[ix, jy, kz] +=
                                wadr * fc^2 * n2r * lr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[ix, jy, kz] += wadr * omir
                    end
                end
            end
        end
    end

    if fc != 0
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            integrals.utheta[ix, jy, kz] =
                thetastrattfc[ix, jy, kz] / fc * integrals.ety[ix, jy, kz]
            integrals.vtheta[ix, jy, kz] =
                -thetastrattfc[ix, jy, kz] / fc * integrals.etx[ix, jy, kz]
        end
    end
end

"""
    compute_gw_integrals!(state::State, wkb_mode::SingleColumn)

Compute gravity wave integrals for single-column WKB mode.

Simplified integration for column models where only vertical momentum
transport matters and horizontal propagation is neglected.

# Arguments

  - `state::State`: Complete simulation state
  - `wkb_mode::SingleColumn`: Single-column mode with vertical-only propagation

# Computed Integrals

Focus on vertical momentum transport:

  - `uw`: Vertical zonal momentum flux `⟨u'w'⟩`
  - `vw`: Vertical meridional momentum flux `⟨v'w'⟩`
  - `etx`, `ety`: Eliassen-Palm flux components
  - `e`: Wave energy density

# Simplified Physics

  - No horizontal momentum fluxes (`uu`, `uv`, `vv` not computed)
  - Only vertical group velocity matters
  - Dispersion relation still includes rotation effects
  - Heat fluxes computed if rotation present

# Efficiency

Much faster than multi-column mode since:

  - No horizontal ray tracking needed
  - Simplified overlap calculations
  - Reduced flux components to compute
"""
function compute_gw_integrals!(state::State, wkb_mode::SingleColumn)
    (; domain, grid) = state
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = domain
    (; dx, dy, dz, x, y, ztildetfc, jac) = grid
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(GWIntegrals)
        getfield(integrals, field) .= 0.0
    end

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + ix] + dx / 2) -
                        max(xr - dxr / 2, x[io + ix] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jy] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jy] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        integrals.uw[ix, jy, kz] +=
                            wadr * kr * cgirz / (1 - (fc / omir)^2)

                        integrals.vw[ix, jy, kz] +=
                            wadr * lr * cgirz / (1 - (fc / omir)^2)

                        if fc != 0
                            integrals.etx[ix, jy, kz] +=
                                wadr * fc^2 * n2r * kr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[ix, jy, kz] +=
                                wadr * fc^2 * n2r * lr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[ix, jy, kz] += wadr * omir
                    end
                end
            end
        end
    end

    if fc != 0
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            integrals.utheta[ix, jy, kz] =
                thetastrattfc[ix, jy, kz] / fc * integrals.ety[ix, jy, kz]
            integrals.vtheta[ix, jy, kz] =
                -thetastrattfc[ix, jy, kz] / fc * integrals.etx[ix, jy, kz]
        end
    end
end

"""
    compute_gw_integrals!(state::State, wkb_mode::SteadyState)

Compute gravity wave integrals for steady-state WKB mode.

Minimal integration for steady-state mountain wave simulations
where only essential momentum flux components are needed.

# Arguments

  - `state::State`: Complete simulation state
  - `wkb_mode::SteadyState`: Steady-state mode with fixed background

# Computed Integrals

Only essential momentum transport:

  - `uw`: Vertical zonal momentum flux
  - `vw`: Vertical meridional momentum flux

# Steady-State Assumptions

  - Background atmosphere is time-independent
  - Wave field reaches quasi-equilibrium
  - Only vertical momentum transport affects mean flow
  - No EP flux computations needed

# Applications

Optimized for:

  - Mountain wave drag parameterizations
  - Long-term climate integrations
  - Cases where transient effects are negligible
"""
function compute_gw_integrals!(state::State, wkb_mode::SteadyState)
    (; domain, grid) = state
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = state.domain
    (; dx, dy, dz, x, y, ztildetfc, jac) = state.grid
    (; sizex, sizey) = state.namelists.domain
    (; branchr) = state.namelists.wkb
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(GWIntegrals)
        getfield(integrals, field) .= 0.0
    end

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min(xr + dxr / 2, x[io + ix] + dx / 2) -
                        max(xr - dxr / 2, x[io + ix] - dx / 2)
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min(yr + dyr / 2, y[jo + jy] + dy / 2) -
                            max(yr - dyr / 2, y[jo + jy] - dy / 2)
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        integrals.uw[ix, jy, kz] += wadr * kr * cgirz

                        integrals.vw[ix, jy, kz] += wadr * lr * cgirz
                    end
                end
            end
        end
    end
end
