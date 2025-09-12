"""
```julia
merge_rays!(state::State)
```

Merge ray volumes by dispatching to a test-case-specific method.

```julia
merge_rays!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
merge_rays!(state::State, testcase::AbstractWKBTestCase)
```

Merge ray volumes by dispatching to a WKB-mode-specific method.

```julia
merge_rays!(state::State, wkb_mode::SteadyState)
```

Return for steady-state WKB mode.

```julia
merge_rays!(state::State, wkb_mode::AbstractWKBMode)
```

Merge ray volumes in grid cells in which their count exceeds a threshold.

This method checks in each grid cell if the number of ray volumes exceeds a maximum that was determined from namelist parameters (`state.wkb.nray_max`). If it does, the ray volumes in that cell are merged such that the new count is smaller or equal to the threshold. This is done by binning them on a spectral grid with logarithmic spacing, defined from the minima and maxima of the contributing negative and positive wavenumbers in all spectral dimensions. The merging is performed such that the bounds of the new ray volumes coincide with the outermost bounds of the old ray volumes and wave action (or wave energy, depending on the merging strategy) is conserved.

# Arguments

  - `state`: Model state.

  - `testcase`: Test case on which the current simulation is based.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.RayOperations.MergedRays`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_spectral_bounds`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_physical_position`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_spectral_position`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_physical_extent`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_spectral_extent`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_surfaces`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_intrinsic_frequency`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_merge_index`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.update_merged_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_wave_action_integral`](@ref)
"""
function merge_rays! end

function merge_rays!(state::State)
    (; testcase) = state.namelists.setting
    merge_rays!(state, testcase)
    return
end

function merge_rays!(state::State, testcase::AbstractTestCase)
    return
end

function merge_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    merge_rays!(state, wkb_mode)
    return
end

function merge_rays!(state::State, wkb_mode::SteadyState)
    return
end

function merge_rays!(state::State, wkb_mode::AbstractWKBMode)
    (; sizex, sizey) = state.namelists.domain
    (; merge_mode) = state.namelists.wkb
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nxray, nyray, nzray, nray_max, nray, rays) = state.wkb

    # Compute ray-volume count before merging.
    @ivy nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    # Initialize array for merged ray volumes.
    merged_rays =
        [MergedRays([[0.0, 0.0] for i in 1:6]..., Ref(0.0)) for i in 1:nray_max]

    # Loop over grid cells.
    @ivy for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if nray[ix, jy, kz] <= nray_max
            continue
        end

        # Set bins in k.
        if sizex > 1
            (kr_min_p, kr_max_p, kr_min_n, kr_max_n) =
                compute_spectral_bounds(rays.k[1:nray[ix, jy, kz], ix, jy, kz])
            dkr_mrg_n = log(kr_max_n / kr_min_n) / (nxray / 2 - 1)
            dkr_mrg_p = log(kr_max_p / kr_min_p) / (nxray / 2 - 1)
        else
            kr_min_p = kr_max_p = kr_min_n = kr_max_n = 0.0
            dkr_mrg_n = dkr_mrg_p = 0.0
        end

        # Set bins in l.
        if sizey > 1
            (lr_min_p, lr_max_p, lr_min_n, lr_max_n) =
                compute_spectral_bounds(rays.l[1:nray[ix, jy, kz], ix, jy, kz])
            dlr_mrg_n = log(lr_max_n / lr_min_n) / (nyray / 2 - 1)
            dlr_mrg_p = log(lr_max_p / lr_min_p) / (nyray / 2 - 1)
        else
            lr_min_p = lr_max_p = lr_min_n = lr_max_n = 0.0
            dlr_mrg_n = dlr_mrg_p = 0.0
        end

        # Set bins in m.
        (mr_min_p, mr_max_p, mr_min_n, mr_max_n) =
            compute_spectral_bounds(rays.m[1:nray[ix, jy, kz], ix, jy, kz])
        dmr_mrg_n = log(mr_max_n / mr_min_n) / (nzray / 2 - 1)
        dmr_mrg_p = log(mr_max_p / mr_min_p) / (nzray / 2 - 1)

        # Reset the merged ray volumes.
        for iray in 1:nray_max
            merged_rays[iray].nr[] = 0.0
        end

        # Loop over ray volumes.
        for iray in 1:nray[ix, jy, kz]
            (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            (dxr, dyr, dzr) = get_physical_extent(rays, (iray, ix, jy, kz))
            (dkr, dlr, dmr) = get_spectral_extent(rays, (iray, ix, jy, kz))
            (axk, ayl, azm) = get_surfaces(rays, (iray, ix, jy, kz))

            nr = rays.dens[iray, ix, jy, kz]
            omegar = compute_intrinsic_frequency(state, (iray, ix, jy, kz))

            # Determine bin index in k.
            if sizex > 1
                fcpspx = axk
                ir_k = compute_merge_index(
                    kr,
                    kr_min_p,
                    kr_max_p,
                    kr_min_n,
                    kr_max_n,
                    dkr_mrg_p,
                    dkr_mrg_n,
                    nxray,
                )
            else
                fcpspx = 1.0
                ir_k = 1
            end

            # Determine bin index in l.
            if sizey > 1
                fcpspy = ayl
                ir_l = compute_merge_index(
                    lr,
                    lr_min_p,
                    lr_max_p,
                    lr_min_n,
                    lr_max_n,
                    dlr_mrg_p,
                    dlr_mrg_n,
                    nyray,
                )
            else
                fcpspy = 1.0
                ir_l = 1
            end

            # Determine bin index in m.
            fcpspz = azm
            ir_m = compute_merge_index(
                mr,
                mr_min_p,
                mr_max_p,
                mr_min_n,
                mr_max_n,
                dmr_mrg_p,
                dmr_mrg_n,
                nzray,
            )

            # Determine flattened bin index.
            if sizex > 1
                if sizey > 1
                    jray =
                        (ir_m - 1) * (nyray - 1) * (nxray - 1) +
                        (ir_l - 1) * (nxray - 1) +
                        ir_k
                else
                    jray = (ir_m - 1) * (nxray - 1) + ir_k
                end
            else
                if sizey > 1
                    jray = (ir_m - 1) * (nyray - 1) + ir_l
                else
                    jray = ir_m
                end
            end

            # Generate the merged ray volumes.
            update_merged_rays!(
                merge_mode,
                merged_rays,
                jray,
                xr,
                dxr,
                yr,
                dyr,
                zr,
                dzr,
                kr,
                dkr,
                lr,
                dlr,
                mr,
                dmr,
                fcpspx,
                fcpspy,
                fcpspz,
                nr,
                omegar,
            )
        end

        # Reset the old ray volumes.
        for field in fieldnames(Rays)
            @. $getfield(rays, field)[:, ix, jy, kz] = 0.0
        end

        # Construct the merged ray volumes.
        iray = 0
        for jray in 1:nray_max
            if merged_rays[jray].nr[] == 0
                continue
            end

            iray += 1

            rays.x[iray, ix, jy, kz] = mean(merged_rays[jray].xr)
            rays.y[iray, ix, jy, kz] = mean(merged_rays[jray].yr)
            rays.z[iray, ix, jy, kz] = mean(merged_rays[jray].zr)

            rays.k[iray, ix, jy, kz] = mean(merged_rays[jray].kr)
            rays.l[iray, ix, jy, kz] = mean(merged_rays[jray].lr)
            rays.m[iray, ix, jy, kz] = mean(merged_rays[jray].mr)

            rays.dxray[iray, ix, jy, kz] = diff(merged_rays[jray].xr)[1]
            rays.dyray[iray, ix, jy, kz] = diff(merged_rays[jray].yr)[1]
            rays.dzray[iray, ix, jy, kz] = diff(merged_rays[jray].zr)[1]

            rays.dkray[iray, ix, jy, kz] = diff(merged_rays[jray].kr)[1]
            rays.dlray[iray, ix, jy, kz] = diff(merged_rays[jray].lr)[1]
            rays.dmray[iray, ix, jy, kz] = diff(merged_rays[jray].mr)[1]

            (axk, ayl, azm) = get_surfaces(rays, (iray, ix, jy, kz))

            omegar = compute_intrinsic_frequency(state, (iray, ix, jy, kz))

            if sizex > 1
                fcpspx = axk
            else
                fcpspx = 1.0
            end

            if sizey > 1
                fcpspy = ayl
            else
                fcpspy = 1.0
            end

            fcpspz = azm

            rays.dens[iray, ix, jy, kz] = compute_wave_action_integral(
                merge_mode,
                merged_rays[jray].nr[],
                1 / omegar,
                1 / fcpspx,
                1 / fcpspy,
                1 / fcpspz,
            )
        end
        nray[ix, jy, kz] = iray
    end

    # Compute ray-volume count after merging.
    @ivy nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_after = MPI.Allreduce(nray_after, +, comm)

    if master && nray_after < nray_before
        println("Number of ray volumes before merging: ", nray_before)
        println("Number of ray volumes after merging: ", nray_after)
        println("")
    end

    return
end
