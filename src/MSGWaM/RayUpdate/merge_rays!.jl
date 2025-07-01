"""
```julia
merge_rays!(state::State)
```

Entry point for ray merging operations based on test case type.

Dispatches to the appropriate merging method depending on the simulation
configuration.

# Arguments

  - `state::State`: Complete simulation state containing ray data
"""
function merge_rays!(state::State)
    (; testcase) = state.namelists.setting
    merge_rays!(state, testcase)
    return
end

"""
```julia
merge_rays!(state::State, testcase::AbstractTestCase)
```

No-op for non-WKB test cases.

Standard test cases don't use ray tracing, so no ray merging is needed.

# Arguments

  - `state::State`: Simulation state (unused)
  - `testcase::AbstractTestCase`: Non-WKB test case
"""
function merge_rays!(state::State, testcase::AbstractTestCase)
    return
end

"""
```julia
merge_rays!(state::State, testcase::AbstractWKBTestCase)
```

Merge rays for WKB test cases based on WKB mode.

Dispatches to the specific WKB mode implementation for ray merging.

# Arguments

  - `state::State`: Simulation state containing WKB configuration and ray data
  - `testcase::AbstractWKBTestCase`: WKB test case specification
"""
function merge_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    merge_rays!(state, wkb_mode)
    return
end

"""
```julia
merge_rays!(state::State, wkb_mode::SteadyState)
```

No-op for steady-state WKB mode.

Steady-state mode typically doesn't require ray merging as the ray
distribution remains relatively stable.

# Arguments

  - `state::State`: Simulation state (unused)
  - `wkb_mode::SteadyState`: Steady-state WKB mode
"""
function merge_rays!(state::State, wkb_mode::SteadyState)
    return
end

"""
```julia
merge_rays!(state::State, wkb_mode::AbstractWKBMode)
```

Merge rays when ray count exceeds maximum per grid cell.

When the number of ray volumes in a grid cell exceeds the maximum allowed,
this function combines nearby rays in phase space to reduce the total count
while preserving wave action and spectral characteristics.

# Arguments

  - `state::State`: Complete simulation state
  - `wkb_mode::AbstractWKBMode`: WKB mode (MultiColumn, SingleColumn, etc.)

# Algorithm

 1. **Check Ray Count**: Only merge if `nray > nray_max` in any cell
 2. **Define Spectral Bins**: Create logarithmic bins in k, l, m space
 3. **Spectral Bounds**: Compute min/max wavenumbers for positive/negative values
 4. **Ray Assignment**: Assign each ray to appropriate spectral bin
 5. **Spatial Merging**: Combine spatial extents (min/max positions)
 6. **Spectral Merging**: Combine spectral extents in each bin
 7. **Wave Action Integration**: Conserve total wave action in each bin
 8. **Reconstruction**: Create merged rays at bin centers

# Spectral Binning

  - **Positive wavenumbers**: Logarithmic spacing from k_min to k_max
  - **Negative wavenumbers**: Logarithmic spacing from -k_max to -k_min
  - **Bin count**: User-specified (nxray, nyray, nzray)

# Wave Action Conservation

Total wave action is preserved: `∑ A_old = ∑ A_new`
where A includes both density and phase space volume factors.

# Benefits

  - Prevents excessive ray counts that would slow computation
  - Maintains spectral representation of wave field
  - Preserves total wave energy and momentum flux
  - Enables long-time integrations

# Statistics

Reports before/after ray counts if merging occurs.
"""
function merge_rays!(state::State, wkb_mode::AbstractWKBMode)
    (; sizex, sizey) = state.namelists.domain
    (; merge_mode) = state.namelists.wkb
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nxray, nyray, nzray, nray_max, nray, rays) = state.wkb

    # Compute ray-volume count before merging.
    @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    # Initialize array for merged ray volumes.
    merged_rays =
        [MergedRays([[0.0, 0.0] for i in 1:6]..., Ref(0.0)) for i in 1:nray_max]

    # Loop over grid cells.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if nray[ix, jy, kz] <= nray_max
            continue
        end

        # Set bins in k.
        if sizex > 1
            @views (kr_min_p, kr_max_p, kr_min_n, kr_max_n) =
                compute_spectral_bounds(rays.k[1:nray[ix, jy, kz], ix, jy, kz])
            dkr_mrg_n = log(kr_max_n / kr_min_n) / (nxray / 2 - 1)
            dkr_mrg_p = log(kr_max_p / kr_min_p) / (nxray / 2 - 1)
        else
            kr_min_p = kr_max_p = kr_min_n = kr_max_n = 0.0
            dkr_mrg_n = dkr_mrg_p = 0.0
        end

        # Set bins in l.
        if sizey > 1
            @views (lr_min_p, lr_max_p, lr_min_n, lr_max_n) =
                compute_spectral_bounds(rays.l[1:nray[ix, jy, kz], ix, jy, kz])
            dlr_mrg_n = log(lr_max_n / lr_min_n) / (nyray / 2 - 1)
            dlr_mrg_p = log(lr_max_p / lr_min_p) / (nyray / 2 - 1)
        else
            lr_min_p = lr_max_p = lr_min_n = lr_max_n = 0.0
            dlr_mrg_n = dlr_mrg_p = 0.0
        end

        # Set bins in m.
        @views (mr_min_p, mr_max_p, mr_min_n, mr_max_n) =
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
            @views getfield(rays, field)[:, ix, jy, kz] .= 0.0
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
    @views nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_after = MPI.Allreduce(nray_after, +, comm)

    if master && nray_after < nray_before
        println("Number of ray volumes before merging: ", nray_before)
        println("Number of ray volumes after merging: ", nray_after)
        println("")
    end

    return
end
