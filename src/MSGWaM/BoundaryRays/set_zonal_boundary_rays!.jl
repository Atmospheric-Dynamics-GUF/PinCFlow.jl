function set_zonal_boundary_rays!(state::State)
    (; namelists, domain) = state
    (; nprocx) = namelists.domain
    (; nxx, sizexx, io, i0, i1, j0, j1, k0, k1) = domain
    (; lx, x) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count.
    set_zonal_boundaries_of_reduced_field!(nray, namelists, domain)

    # Set ray-volumes properties.
    if nprocx > 1
        set_zonal_halo_rays!(state)
    else
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1)
            if nray[i0 - 1, jy, kz] > 0
                for iray in 1:nray[i0 - 1, jy, kz]
                    copy_rays!(rays, (iray, i1, jy, kz), (iray, i0 - 1, jy, kz))
                end
            end

            if nray[i1 + 1, jy, kz] > 0
                for iray in 1:nray[i1 + 1, jy, kz]
                    copy_rays!(rays, (iray, i0, jy, kz), (iray, i1 + 1, jy, kz))
                end
            end
        end
    end

    if io == 0
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1)
            for ix in (i0 - 1):i0
                if nray[ix, jy, kz] > 0
                    for iray in 1:nray[ix, jy, kz]
                        xr = rays.x[iray, ix, jy, kz]
                        xrt = xr - lx[2] + lx[1]

                        if abs(xrt - x[io + ix]) < abs(xr - x[io + ix])
                            xr = xrt
                        end

                        rays.x[iray, ix, jy, kz] = xr
                    end
                end
            end
        end
    end

    if io + nxx == sizexx
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1)
            for ix in i1:(i1 + 1)
                if nray[ix, jy, kz] > 0
                    for iray in 1:nray[ix, jy, kz]
                        xr = rays.x[iray, ix, jy, kz]
                        xrt = xr + lx[2] - lx[1]

                        if abs(xrt - x[io + ix]) < abs(xr - x[io + ix])
                            xr = xrt
                        end

                        rays.x[iray, ix, jy, kz] = xr
                    end
                end
            end
        end
    end

    return
end
