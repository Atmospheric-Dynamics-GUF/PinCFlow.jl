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

function merge_rays!(state::State, wkb_mode::Union{SingleColumn, MultiColumn})
    (; sizex, sizey) = state.namelists.domain
    (; merge_mode) = state.namelists.wkb
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nxray, nyray, nzray, nray_max, nray, rays) = state.wkb

    # Compute ray-volume count before merging.
    @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    # Initialize array for merged ray volumes.
    merged_rays = [MergedRays() for i in 1:nray_max]

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
            merged_rays[jray] = MergedRays(
                merge_mode,
                merged_rays[jray],
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

        # Replace ray volumes.
        iray = 0
        for jray in 1:nray_max
            if merged_rays[jray].nr == 0
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

            (axk, ayl, azm) =
                get_surfaces(rays, (iray, ix, jy, kz))

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

            if merge_mode == ConstantWaveAction()
                rays.dens[iray, ix, jy, kz] =
                    merged_rays[jray].nr / (fcpspx * fcpspy * fcpspz)
            elseif merge_mode == ConstantWaveEnergy()
                rays.dens[iray, ix, jy, kz] =
                    merged_rays[jray].nr / (omegar * fcpspx * fcpspy * fcpspz)
            end
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
