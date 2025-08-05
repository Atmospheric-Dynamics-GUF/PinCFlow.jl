"""
```julia
initialize_rays!(state::State)
```

Complete the initialization of MSGWaM by dispatching to a test-case-specific method.

```julia
initialize_rays!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
initialize_rays!(state::State, testcase::AbstractWKBTestCase)
```

Complete the initialization of MSGWaM for WKB test cases.

In each grid cell, `nwm` wave modes are computed, using e.g. `activate_orographic_source!` for mountain waves. For each of these modes, `nrxl * nryl * nrzl * nrk * nrl * nrm` ray volumes are then defined such that they evenly divide the volume one would get for `nrxl = nryl = nrzl = nrk = nrl = nrm = 1` (the parameters are taken from `state.namelists.wkb`). Finally, the maximum group velocities are determined for the corresponding CFL condition that is used in the computation of the time step.

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)
"""
function initialize_rays! end

function initialize_rays!(state::State)
    (; testcase) = state.namelists.setting
    initialize_rays!(state, testcase)
    return
end

function initialize_rays!(state::State, testcase::AbstractTestCase)
    return
end

function initialize_rays!(state::State, testcase::AbstractWKBTestCase)
    (; sizex, sizey, sizez) = state.namelists.domain
    (; testcase) = state.namelists.setting
    (; coriolis_frequency) = state.namelists.atmosphere
    (;
        xrmin_dim,
        xrmax_dim,
        yrmin_dim,
        yrmax_dim,
        nrxl,
        nryl,
        nrzl,
        nrk_init,
        nrl_init,
        nrm_init,
        nwm,
        fac_dk_init,
        fac_dl_init,
        fac_dm_init,
        wkb_mode,
        nwm,
    ) = state.namelists.wkb
    (; lref, tref) = state.constants
    (; comm, master, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; lx, ly, lz, dx, dy, dz, x, y, ztfc, jac) = state.grid
    (;
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    # Non-dimensionalize boundaries for ray-volume propagation.
    xrmin = xrmin_dim / lref
    xrmax = xrmax_dim / lref
    yrmin = yrmin_dim / lref
    yrmax = yrmax_dim / lref

    # Set zonal index bounds.
    if testcase == WKBMountainWave()
        ix0 = i0
        ix1 = i1
    else
        ix0 = max(i0, floor(Int, (xrmin - lx[1]) / dx) + i0 - io)
        ix1 = min(i1, floor(Int, (xrmax - lx[1]) / dx) + i0 - io)
    end

    # Set meridional index bounds.
    if testcase == WKBMountainWave()
        jy0 = j0
        jy1 = j1
    else
        jy0 = max(j0, floor(Int, (yrmin - ly[1]) / dy) + j0 - jo)
        jy1 = min(j1, floor(Int, (yrmax - ly[1]) / dy) + j0 - jo)
    end

    # Set vertical index bounds.
    if testcase == WKBMountainWave() && ko == 0
        kz0 = k0 - 1
        kz1 = k0 - 1
    else
        kz0 = k0
        kz1 = k1
    end

    # Initialize local arrays.
    omi_ini = zeros(nwm, nxx, nyy, nzz)
    wnk_ini = zeros(nwm, nxx, nyy, nzz)
    wnl_ini = zeros(nwm, nxx, nyy, nzz)
    wnm_ini = zeros(nwm, nxx, nyy, nzz)
    wad_ini = zeros(nwm, nxx, nyy, nzz)

    if wkb_mode == SteadyState() && testcase != WKBMountainWave()
        error(
            "Error in initialize_rays!: SteadyState is implemented for WKBMountainWave only!",
        )
    end

    if testcase == WKBMountainWave()
        activate_orographic_source!(
            state,
            omi_ini,
            wnk_ini,
            wnl_ini,
            wnm_ini,
            wad_ini,
        )
    end

    dk_ini_nd = 0.0
    dl_ini_nd = 0.0
    dm_ini_nd = 0.0

    # Loop over all grid cells with ray volumes.
    for kz in kz0:kz1, jy in jy0:jy1, ix in ix0:ix1
        iray = 0
        i_sfc = 0

        # Loop over all ray volumes within a spatial cell.
        for ix2 in 1:nrxl,
            ik in 1:nrk_init,
            jy2 in 1:nryl,
            jl in 1:nrl_init,
            kz2 in 1:nrzl,
            km in 1:nrm_init,
            iwm in 1:nwm

            # Set ray-volume indices.
            if testcase == WKBMountainWave()
                i_sfc += 1

                # Set surface indices.
                surface_indices.ix2_sfc[i_sfc] = ix2
                surface_indices.jy2_sfc[i_sfc] = jy2
                surface_indices.kz2_sfc[i_sfc] = kz2
                surface_indices.ik_sfc[i_sfc] = ik
                surface_indices.jl_sfc[i_sfc] = jl
                surface_indices.km_sfc[i_sfc] = km
                surface_indices.iwm_sfc[i_sfc] = iwm

                # Set surface ray-volume index.
                if wad_ini[iwm, ix, jy, kz] == 0.0
                    surface_indices.ir_sfc[i_sfc, ix, jy] = -1
                    continue
                else
                    iray += 1
                    surface_indices.ir_sfc[i_sfc, ix, jy] = iray
                end
            else
                iray += 1
            end

            # Set ray-volume positions.
            rays.x[iray, ix, jy, kz] =
                (x[io + ix] - 0.5 * dx + (ix2 - 0.5) * dx / nrxl)
            rays.y[iray, ix, jy, kz] =
                (y[jo + jy] - 0.5 * dy + (jy2 - 0.5) * dy / nryl)
            rays.z[iray, ix, jy, kz] = (
                ztfc[ix, jy, kz] - 0.5 * jac[ix, jy, kz] * dz +
                (kz2 - 0.5) * jac[ix, jy, kz] * dz / nrzl
            )

            xr = rays.x[iray, ix, jy, kz]
            yr = rays.y[iray, ix, jy, kz]
            zr = rays.z[iray, ix, jy, kz]

            # Check if ray volume is too low.
            if zr < lz[1] - dz
                error(
                    "Error in initialize_rays!: Ray volume",
                    iray,
                    "at",
                    ix,
                    jy,
                    kz,
                    "is too low!",
                )
            end

            # Compute local stratification.
            n2r = interpolate_stratification(zr, state, N2())

            # Set spatial extents.
            rays.dxray[iray, ix, jy, kz] = dx / nrxl
            rays.dyray[iray, ix, jy, kz] = dy / nryl
            rays.dzray[iray, ix, jy, kz] = jac[ix, jy, kz] * dz / nrzl

            wnk0 = wnk_ini[iwm, ix, jy, kz]
            wnl0 = wnl_ini[iwm, ix, jy, kz]
            wnm0 = wnm_ini[iwm, ix, jy, kz]

            # Ensure correct wavenumber extents.
            if testcase == WKBMountainWave() && sizex > 1
                dk_ini_nd = fac_dk_init * sqrt(wnk0^2 + wnl0^2)
            end
            if testcase == WKBMountainWave() && sizey > 1
                dl_ini_nd = fac_dl_init * sqrt(wnk0^2 + wnl0^2)
            end
            if wnm0 == 0.0
                error("Error in WKB: wnm0 = 0!")
            else
                dm_ini_nd = fac_dm_init * abs(wnm0)
            end

            # Set ray-volume wavenumbers.
            rays.k[iray, ix, jy, kz] =
                (wnk0 - 0.5 * dk_ini_nd + (ik - 0.5) * dk_ini_nd / nrk_init)
            rays.l[iray, ix, jy, kz] =
                (wnl0 - 0.5 * dl_ini_nd + (jl - 0.5) * dl_ini_nd / nrl_init)
            rays.m[iray, ix, jy, kz] =
                (wnm0 - 0.5 * dm_ini_nd + (km - 0.5) * dm_ini_nd / nrm_init)

            # Set spectral extents.
            rays.dkray[iray, ix, jy, kz] = dk_ini_nd / nrk_init
            rays.dlray[iray, ix, jy, kz] = dl_ini_nd / nrl_init
            rays.dmray[iray, ix, jy, kz] = dm_ini_nd / nrm_init

            # Set spectral volume.
            pspvol = dm_ini_nd
            if sizex > 1
                pspvol = pspvol * dk_ini_nd
            end
            if sizey > 1
                pspvol = pspvol * dl_ini_nd
            end

            # Set phase-space wave-action density.
            if kz == sizez
                rays.dens[iray, ix, jy, kz] = 0.0
            else
                rays.dens[iray, ix, jy, kz] = wad_ini[iwm, ix, jy, kz] / pspvol
            end

            # Interpolate winds to ray-volume position.
            uxr = interpolate_mean_flow(xr, yr, zr, state, U())
            vyr = interpolate_mean_flow(xr, yr, zr, state, V())
            wzr = interpolate_mean_flow(xr, yr, zr, state, W())

            wnrk = rays.k[iray, ix, jy, kz]
            wnrl = rays.l[iray, ix, jy, kz]
            wnrm = rays.m[iray, ix, jy, kz]
            wnrh = sqrt(wnrk^2 + wnrl^2)
            omir = omi_ini[iwm, ix, jy, kz]

            # Compute maximum group velocities.
            cgirx = wnrk * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(uxr + cgirx) > abs(cgx_max[])
                cgx_max[] = abs(uxr + cgirx)
            end
            cgiry = wnrl * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(vyr + cgiry) > abs(cgy_max[])
                cgy_max[] = abs(vyr + cgiry)
            end
            cgirz = -wnrm * (omir^2 - fc^2) / (omir * (wnrh^2 + wnrm^2))
            if abs(wzr + cgirz) > abs(cgz_max[ix, jy, kz])
                cgz_max[ix, jy, kz] = max(cgz_max[ix, jy, kz], abs(wzr + cgirz))
            end
        end

        # Set ray-volume count.
        nray[ix, jy, kz] = iray
        if iray > nray_wrk
            error(
                "Error in initialize_rays!: nray",
                [ix, jy, kz],
                " > nray_wrk =",
                nray_wrk,
            )
        end

        # Check if surface ray-volume count is correct.
        if testcase == WKBMountainWave()
            if i_sfc != n_sfc
                error(
                    "Error in initialize_rays!: i_sfc =",
                    i_sfc,
                    "/= n_sfc =",
                    n_sfc,
                    "at (ix, jy, kz) = ",
                    (ix, jy, kz),
                )
            end
        end
    end

    # Compute global ray-volume count.
    local_sum = sum(nray[ix0:ix1, jy0:jy1, kz0:kz1])
    global_sum = MPI.Allreduce(local_sum, +, comm)

    # Print information.
    if master
        println("MSGWaM:")
        println("Global ray-volume count: ", global_sum)
        println("Maximum number of ray volumes per cell: ", nray_max)
        println("")
    end

    return
end
