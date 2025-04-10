function set_meridional_boundary_ray_volumes!(state::State)
    (; namelists, domain) = state
    (; sizey, nprocy) = namelists.domain
    (; nyy, jo, i0, i1, j0, j1, k0, k1) = domain
    (; ly, y) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count.
    set_zonal_boundaries_of_reduced_field!(nray, namelists, domain)

    # Set ray-volumes properties.
    if nprocy > 1
        error("Halo routines are not working yet!")
        set_zonal_halo_ray_volumes!(state)
    else
        for kz in (k0 - 1):(k1 + 1), ix in (i0 - 1):(i1 + 1)
            if nray[ix, j0 - 1, kz] > 0
                for iray in 1:nray[ix, j0 - 1, kz]
                    copy_ray_volume!(
                        rays,
                        (iray, ix, j1, kz),
                        (iray, ix, j0 - 1, kz),
                    )
                end
            end

            if nray[ix, j1 + 1, kz] > 0
                for iray in 1:nray[ix, j1 + 1, kz]
                    copy_ray_volume!(
                        rays,
                        (iray, ix, j0, kz),
                        (iray, ix, j1 + 1, kz),
                    )
                end
            end
        end
    end

    if jo + nyy == sizey
        for kz in (k0 - 1):(k1 + 1), ix in (i0 - 1):(i1 + 1)
            for jy in (j0 - 1):j0
                if nray[ix, jy, kz] > 0
                    for iray in 1:nray[ix, jy, kz]
                        yr = rays.y[iray, ix, jy, kz]
                        yrt = yr - ly[2] + ly[1]

                        if abs(yrt - y[jo + jy]) < abs(yr - y[jo + jy])
                            yr = yrt
                        end

                        rays.y[iray, ix, jy, kz] = yr
                    end
                end
            end
        end
    end

    if io == 0
        for kz in (k0 - 1):(k1 + 1), ix in (i0 - 1):(i1 + 1)
            for jy in j1:(j1 + 1)
                if nray[ix, jy, kz] > 0
                    for iray in 1:nray[ix, jy, kz]
                        yr = rays.y[iray, ix, jy, kz]
                        yrt = yr + ly[2] - ly[1]

                        if abs(yrt - y[jo + jy]) < abs(yr - y[jo + jy])
                            yr = yrt
                        end

                        rays.y[iray, ix, jy, kz] = yr
                    end
                end
            end
        end
    end

    return
end
