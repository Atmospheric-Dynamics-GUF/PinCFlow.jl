"""
```julia
merge_rays!(state::State)
```

Merge ray volumes by dispatching to a WKB-mode-specific method.

```julia
merge_rays!(state::State, wkb_mode::Union{NoWKB, SteadyState})
```

Return for configurations without WKB / with steady-state WKB.

```julia
merge_rays!(state::State, wkb_mode::Union{SingleColumn, MultiColumn})
```

Merge ray volumes in grid cells in which their count exceeds a threshold.

This method checks in each grid cell if the number of ray volumes exceeds a maximum that was determined from namelist parameters (`state.wkb.nray_max`). If it does, the ray volumes in that cell are merged such that the new count is smaller or equal to the threshold. This is done by binning them on a spectral grid with logarithmic spacing, defined from the minima and maxima of the contributing negative and positive wavenumbers in all spectral dimensions. The merging is performed such that the bounds of the new ray volumes coincide with the outermost bounds of the old ray volumes and wave action (or wave energy, depending on the merging strategy) is conserved.

# Arguments

  - `state`: Model state.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

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
    (; wkb_mode) = state.namelists.wkb
    merge_rays!(state, wkb_mode)
    return
end

function merge_rays!(state::State, wkb_mode::Union{NoWKB, SteadyState})
    return
end

function merge_rays!(state::State, wkb_mode::Union{SingleColumn, MultiColumn})
    (; x_size, y_size) = state.namelists.domain
    (; merge_mode) = state.namelists.wkb
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nxray, nyray, nzray, nray_max, nray, rays, merged_rays) = state.wkb

    # Compute ray-volume count before merging.
    @ivy nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    # Loop over grid cells.
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        if nray[i, j, k] <= nray_max
            continue
        end

        # Set bins in k.
        if x_size > 1
            (kr_min_p, kr_max_p, kr_min_n, kr_max_n) =
                compute_spectral_bounds(rays.k[1:nray[i, j, k], i, j, k])
            dkr_mrg_n = log(kr_max_n / kr_min_n) / (nxray / 2 - 1)
            dkr_mrg_p = log(kr_max_p / kr_min_p) / (nxray / 2 - 1)
        else
            kr_min_p = kr_max_p = kr_min_n = kr_max_n = 0.0
            dkr_mrg_n = dkr_mrg_p = 0.0
        end

        # Set bins in l.
        if y_size > 1
            (lr_min_p, lr_max_p, lr_min_n, lr_max_n) =
                compute_spectral_bounds(rays.l[1:nray[i, j, k], i, j, k])
            dlr_mrg_n = log(lr_max_n / lr_min_n) / (nyray / 2 - 1)
            dlr_mrg_p = log(lr_max_p / lr_min_p) / (nyray / 2 - 1)
        else
            lr_min_p = lr_max_p = lr_min_n = lr_max_n = 0.0
            dlr_mrg_n = dlr_mrg_p = 0.0
        end

        # Set bins in m.
        (mr_min_p, mr_max_p, mr_min_n, mr_max_n) =
            compute_spectral_bounds(rays.m[1:nray[i, j, k], i, j, k])
        dmr_mrg_n = log(mr_max_n / mr_min_n) / (nzray / 2 - 1)
        dmr_mrg_p = log(mr_max_p / mr_min_p) / (nzray / 2 - 1)

        # Reset the merged ray volumes.
        merged_rays.nr .= 0.0

        # Loop over ray volumes.
        for r in 1:nray[i, j, k]
            (xr, yr, zr) = get_physical_position(rays, r, i, j, k)
            (kr, lr, mr) = get_spectral_position(rays, r, i, j, k)
            (dxr, dyr, dzr) = get_physical_extent(rays, r, i, j, k)
            (dkr, dlr, dmr) = get_spectral_extent(rays, r, i, j, k)
            (axk, ayl, azm) = get_surfaces(rays, r, i, j, k)

            nr = rays.dens[r, i, j, k]
            omegar = compute_intrinsic_frequency(state, r, i, j, k)

            # Determine bin index in k.
            if x_size > 1
                fcpspx = axk
                rk = compute_merge_index(
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
                rk = 1
            end

            # Determine bin index in l.
            if y_size > 1
                fcpspy = ayl
                rl = compute_merge_index(
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
                rl = 1
            end

            # Determine bin index in m.
            fcpspz = azm
            rm = compute_merge_index(
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
            if x_size > 1
                if y_size > 1
                    bin =
                        (rm - 1) * (nyray - 1) * (nxray - 1) +
                        (rl - 1) * (nxray - 1) +
                        rk
                else
                    bin = (rm - 1) * (nxray - 1) + rk
                end
            else
                if y_size > 1
                    bin = (rm - 1) * (nyray - 1) + rl
                else
                    bin = rm
                end
            end

            # Generate the merged ray volumes.
            update_merged_rays!(
                merge_mode,
                merged_rays,
                bin,
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
            getfield(rays, field)[:, i, j, k] .= 0.0
        end

        # Construct the merged ray volumes.
        r = 0
        for bin in 1:nray_max
            if merged_rays.nr[bin] == 0
                continue
            end

            r += 1

            rays.x[r, i, j, k] = sum(merged_rays.xr[:, bin]) / 2
            rays.y[r, i, j, k] = sum(merged_rays.yr[:, bin]) / 2
            rays.z[r, i, j, k] = sum(merged_rays.zr[:, bin]) / 2

            rays.k[r, i, j, k] = sum(merged_rays.kr[:, bin]) / 2
            rays.l[r, i, j, k] = sum(merged_rays.lr[:, bin]) / 2
            rays.m[r, i, j, k] = sum(merged_rays.mr[:, bin]) / 2

            rays.dxray[r, i, j, k] = diff(merged_rays.xr[:, bin])[1]
            rays.dyray[r, i, j, k] = diff(merged_rays.yr[:, bin])[1]
            rays.dzray[r, i, j, k] = diff(merged_rays.zr[:, bin])[1]

            rays.dkray[r, i, j, k] = diff(merged_rays.kr[:, bin])[1]
            rays.dlray[r, i, j, k] = diff(merged_rays.lr[:, bin])[1]
            rays.dmray[r, i, j, k] = diff(merged_rays.mr[:, bin])[1]

            #CHANGES
            # set phase to zero after merging
            rays.dphi[iray, i, j, k] = 0.0

            (axk, ayl, azm) = get_surfaces(rays, r, i, j, k)

            omegar = compute_intrinsic_frequency(state, r, i, j, k)

            if x_size > 1
                fcpspx = axk
            else
                fcpspx = 1.0
            end

            if y_size > 1
                fcpspy = ayl
            else
                fcpspy = 1.0
            end

            fcpspz = azm

            rays.dens[r, i, j, k] = compute_wave_action_integral(
                merge_mode,
                merged_rays.nr[bin],
                1 / omegar,
                1 / fcpspx,
                1 / fcpspy,
                1 / fcpspz,
            )
        end
        nray[i, j, k] = r
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
