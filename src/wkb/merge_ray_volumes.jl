function merge_ray_volumes!(state::State)
    (; wkb_mode) = state.namelists.wkb
    merge_ray_volumes(state, wkb_mode)
    return
end

function merge_ray_volumes!(state::State, wkb_mode::SteadyState)
    return
end

function merge_ray_volumes!(
    state::State,
    wkb_mode::Union{SingleColumn, MultiColumn},
)

    # Compute ray-volume count before merging.
    @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    # Initialize array for merged ray volumes.
    merged_rays = [MergedRayVolume() for i in 1:nray_max]

    # Loop over grid cells.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if nray[ix, jy, kz] <= nray_max
            continue
        end

        # Set bins in k.
        if sizex > 1
            @views (wnrk_min_p, wnrk_max_p, wnrk_min_n, wnrk_max_n) =
                compute_spectral_bounds(rays.k[1:nray[ix, jy, kz], ix, jy, kz])
            dwnrk_mrg_n = log(wnrk_max_n / wnrk_min_n) / (nxray / 2 - 1)
            dwnrk_mrg_p = log(wnrk_max_p / wnrk_min_p) / (nxray / 2 - 1)
        else
            wnrk_min_p = wnrk_max_p = wnrk_min_n = wnrk_max_n = 0.0
            dwnrk_mrg_n = dwnrk_mrg_p = 0.0
        end

        # Set bins in l.
        if sizey > 1
            @views (wnrl_min_p, wnrl_max_p, wnrl_min_n, wnrl_max_n) =
                compute_spectral_bounds(rays.l[1:nray[ix, jy, kz], ix, jy, kz])
            dwnrl_mrg_n = log(wnrl_max_n / wnrl_min_n) / (nyray / 2 - 1)
            dwnrl_mrg_p = log(wnrl_max_p / wnrl_min_p) / (nyray / 2 - 1)
        else
            wnrl_min_p = wnrl_max_p = wnrl_min_n = wnrl_max_n = 0.0
            dwnrl_mrg_n = dwnrl_mrg_p = 0.0
        end

        # Set bins in m.
        @views (wnrm_min_p, wnrm_max_p, wnrm_min_n, wnrm_max_n) =
            compute_spectral_bounds(rays.m[1:nray[ix, jy, kz], ix, jy, kz])
        dwnrm_mrg_n = log(wnrm_max_n / wnrm_min_n) / (nzray / 2 - 1)
        dwnrm_mrg_p = log(wnrm_max_p / wnrm_min_p) / (nzray / 2 - 1)

        # Loop over ray volumes.
        for iray in 1:nray[ix, jy, kz]
            xr = rays.x[iray, ix, jy, kz]
            dxr = rays.dxray[iray, ix, jy, kz]

            yr = rays.y[iray, ix, jy, kz]
            dyr = rays.dyray[iray, ix, jy, kz]

            zr = rays.z[iray, ix, jy, kz]
            dzr = rays.dzray[iray, ix, jy, kz]

            wnrk = rays.k[iray, ix, jy, kz]
            dwnrk = rays.dkray[iray, ix, jy, kz]

            wnrl = rays.l[iray, ix, jy, kz]
            dwnrl = rays.dlray[iray, ix, jy, kz]

            wnrm = rays.m[iray, ix, jy, kz]
            dwnrm = rays.dmray[iray, ix, jy, kz]

            axk = rays.area_xk[iray, ix, jy, kz]
            ayl = rays.area_yl[iray, ix, jy, kz]
            azm = rays.area_zm[iray, ix, jy, kz]

            wdr = rays.dens[iray, ix, jy, kz]
            omir = ray.omega[iray, ix, jy, kz]

            # Determine bin index in k.
            if sizex > 1
                fcpspx = axk
                ir_k = compute_merge_index(
                    wnrk,
                    wnrk_min_p,
                    wnrk_max_p,
                    wnrk_min_n,
                    wnrk_max_n,
                    dwnrk_mrg_p,
                    dwnrk_mrg_n,
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
                    wnrl,
                    wnrl_min_p,
                    wnrl_max_p,
                    wnrl_min_n,
                    wnrl_max_n,
                    dwnrl_mrg_p,
                    dwnrl_mrg_n,
                    nyray,
                )
            else
                fcpspy = 1.0
                ir_l = 1
            end

            # Determine bin index in m.
            fcpspz = azm
            ir_m = compute_merge_index(
                wnrm,
                wnrm_min_p,
                wnrm_max_p,
                wnrm_min_n,
                wnrm_max_n,
                dwnrm_mrg_p,
                dwnrm_mrg_n,
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
            merged_rays[jray] = MergedRayVolume(
                merge_mode,
                merged_rays[jray],
                xr,
                dxr,
                yr,
                dyr,
                zr,
                dzr,
                wnrk,
                dwnrk,
                wnrl,
                dwnrl,
                wnrm,
                dwnrm,
                fcpspx,
                fcpspy,
                fcpspz,
                wdr,
                omir,
            )
        end

        # Replace ray volumes.
        iray = 0
        for jray in 1:nray_max
            if merged_rays[jray].dens == 0
                continue
            end

            iray += 1

            rays.x[iray, ix, jy, kz] =
                (merged_rays[jray].xmax + merged_rays[jray].xmin) / 2
            rays.dxray[iray, ix, jy, kz] =
                merged_rays[jray].xmax - merged_rays[jray].xmin

            rays.y[iray, ix, jy, kz] =
                (merged_rays[jray].ymax + merged_rays[jray].ymin) / 2
            rays.dyray[iray, ix, jy, kz] =
                merged_rays[jray].ymax - merged_rays[jray].ymin

            rays.z[iray, ix, jy, kz] =
                (merged_rays[jray].zmax + merged_rays[jray].zmin) / 2
            rays.dzray[iray, ix, jy, kz] =
                merged_rays[jray].zmax - merged_rays[jray].zmin

            rays.k[iray, ix, jy, kz] =
                (merged_rays[jray].kmax + merged_rays[jray].kmin) / 2
            rays.dkray[iray, ix, jy, kz] =
                merged_rays[jray].kmax - merged_rays[jray].kmin

            rays.l[iray, ix, jy, kz] =
                (merged_rays[jray].lmax + merged_rays[jray].lmin) / 2
            rays.dlray[iray, ix, jy, kz] =
                merged_rays[jray].lmax - merged_rays[jray].l_min_n

            rays.m[iray, ix, jy, kz] =
                (merged_rays[jray].mmax + merged_rays[jray].mmin) / 2
            rays.dmray[iray, ix, jy, kz] =
                merged_rays[jray].mmax - merged_rays[jray].mmin

            rays.omega[iray, ix, jy, kz] =
                compute_intrinsic_frequency(state, (iray, ix, jy, kz))

            (
                rays.area_xk[iray, ix, jy, kz],
                rays.area_yl[iray, ix, jy, kz],
                rays.area_zm[iray, ix, jy, kz],
            ) = compute_ray_volumes_surfaces(rays, (iray, ix, jy, kz))

            if sizex > 1
                fcpspx = rays.area_xk[iray, ix, jy, kz]
            else
                fcpspx = 1.0
            end

            if sizey > 1
                fcpspy = rays.area_yl[iray, ix, jy, kz]
            else
                fcpspy = 1.0
            end

            fcpspz = rays.area_zm[iray, ix, jy, kz]

            if merge_mode == ConstantWaveAction()
                rays.dens[iray, ix, jy, kz] =
                    merged_rays[jray].dens / (fcpspx * fcpspy * fcpspz)
            elseif merge_mode == ConstantWaveEnergy()
                rays.dens[iray, ix, jy, kz] =
                    merged_rays[jray].dens / (omir * fcpspx * fcpspy * fcpspz)
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
    end

    return
end