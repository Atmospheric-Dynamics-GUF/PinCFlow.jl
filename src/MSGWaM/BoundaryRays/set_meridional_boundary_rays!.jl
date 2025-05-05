function set_meridional_boundary_rays!(state::State)
    (; namelists, domain) = state
    (; npy) = namelists.domain
    (; nyy, sizeyy, jo, i0, i1, j0, j1, k0, k1) = domain
    (; ly, y) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count.
    set_meridional_boundaries_of_field!(
        nray,
        namelists,
        domain;
        layers = (1, 1, 1),
    )

    # Set ray-volumes properties.
    if npy > 1
        set_meridional_halo_rays!(state)
    else
        for kz in (k0 - 1):(k1 + 1), ix in (i0 - 1):(i1 + 1)
            for iray in 1:nray[ix, j0 - 1, kz]
                copy_rays!(rays, (iray, ix, j1, kz), (iray, ix, j0 - 1, kz))
            end

            for iray in 1:nray[ix, j1 + 1, kz]
                copy_rays!(rays, (iray, ix, j0, kz), (iray, ix, j1 + 1, kz))
            end
        end
    end

    if jo == 0
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):j0, ix in (i0 - 1):(i1 + 1)
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

    if jo + nyy == sizeyy
        for kz in (k0 - 1):(k1 + 1), jy in j1:(j1 + 1), ix in (i0 - 1):(i1 + 1)
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

    return
end
