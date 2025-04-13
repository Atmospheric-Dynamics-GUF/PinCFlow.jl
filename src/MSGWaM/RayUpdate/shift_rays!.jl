function shift_rays!(state::State)
    (; wkb_mode) = state.namelists.wkb
    shift_rays!(state, wkb_mode)
    return
end

function shift_rays!(state::State, wkb_mode::SteadyState)
    return
end

function shift_rays!(state::State, wkb_mode::SingleColumn)
    (; zboundaries) = state.namelists.setting

    set_vertical_boundary_rays!(state, zboundaries)
    shift_rays!(state, Z())
    set_vertical_boundary_rays!(state, zboundaries)
    remove_rays!(state)

    check_ray_volumes(state)

    return
end

function shift_rays!(state::State, wkb_mode::MultiColumn)
    (; sizex, sizey) = state.namelists.domain
    (; zboundaries) = state.namelists.setting

    if sizex > 1
        set_zonal_boundary_rays!(state)
        shift_rays!(state, X())
        set_zonal_boundary_rays!(state)
        remove_rays!(state)
    end

    if sizey > 1
        set_meridional_boundary_rays!(state)
        shift_rays!(state, Y())
        set_meridional_boundary_rays!(state)
        remove_rays!(state)
    end

    set_vertical_boundary_rays!(state, zboundaries)
    shift_rays!(state, Z())
    set_vertical_boundary_rays!(state, zboundaries)
    remove_rays!(state)

    check_rays(state)

    return
end

function shift_rays!(state::State, direction::X)
    (; io, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, dx) = state.grid
    (; nray, rays) = state.wkb

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] != 0.0
                xr = rays.x[iray, ixrv, jyrv, kzrv]
                ix = round(Int, (xr - lx[1] - dx / 2) / dx) + i0 - io

                if ix != ixrv
                    if i0 <= ix <= i1
                        nray[ix, jyrv, kzrv] += 1
                        jray = nray[ix, jyrv, kzrv]
                        if jray > nray_wrk
                            error("Error in shift_rays!: nray > nray_wrk!")
                        end
                        copy_ray_volume!(
                            rays,
                            (iray, ixrv, jyrv, kzrv),
                            (jray, ix, jyrv, kzrv),
                        )
                    end
                    rays.dens[iray, ix, jyrv, kzrv] = 0.0
                end
            end
        end
    end
end

function shift_rays!(state::State, direction::Y)
    (; jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; ly, dy) = state.grid
    (; nray, rays) = state.wkb

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] != 0.0
                yr = rays.y[iray, ixrv, jyrv, kzrv]
                jy = round(Int, (yr - ly[1] - dy / 2) / dy) + j0 - jo

                if jy != jyrv
                    if j0 <= jy <= j1
                        nray[ixrv, jy, kzrv] += 1
                        jray = nray[ixrv, jy, kzrv]
                        if jray > nray_wrk
                            error("Error in shift_rays!: nray > nray_wrk!")
                        end
                        copy_ray_volume!(
                            rays,
                            (iray, ixrv, jyrv, kzrv),
                            (jray, ixrv, jy, kzrv),
                        )
                    end
                    rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
                end
            end
        end
    end
end

function shift_rays!(state::State, direction::Z)
    (; domain, grid) = state
    (; i0, i1, j0, j1, k0, k1) = domain
    (; nray, rays) = state.wkb

    for kzrv in k0:k1, jyrv in (j0 - 1):(j1 + 1), ixrv in (i0 - 1):(i1 + 1)
        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] != 0.0
                zr = rays.z[iray, ixrv, jyrv, kzrv]
                kz = get_next_half_level(ixrv, jyrv, zr, domain, grid)

                if kz > k1
                    kz = k1
                end
                if kz < k0
                    kz = k0
                end

                if kz != kzrv
                    nray[ixrv, jyrv, kz] += 1
                    jray = nray[ixrv, jyrv, kz]
                    if jray > nray_wrk
                        error("Error in shift_rays!: nray > nray_wrk!")
                    end
                    copy_rays!(
                        rays,
                        (iray, ixrv, jyrv, kzrv),
                        (jray, ixrv, jyrv, kz),
                    )
                    rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
                end
            end
        end
    end

    return
end
