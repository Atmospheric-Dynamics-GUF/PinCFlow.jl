function shift_rays!(state::State)
    (; testcase) = state.namelists.setting
    shift_rays!(state, testcase)
    return
end

function shift_rays!(state::State, testcase::AbstractTestCase)
    return
end

function shift_rays!(state::State, testcase::AbstractWKBTestCase)
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

    check_rays(state)

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
    (; sizezz, nzz, io, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, dx) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for kzrv in kz0:kz1, jyrv in (j0 - 1):(j1 + 1), ixrv in (i0 - 1):(i1 + 1)
        for iray in 1:nray[ixrv, jyrv, kzrv]
            xr = rays.x[iray, ixrv, jyrv, kzrv]
            ix = floor(Int, (xr - lx[1]) / dx) + i0 - io

            if ix != ixrv
                if abs(ix - ixrv) > 1
                    error("Error in shift_rays!: abs(ix - ixrv) > 1!")
                end
                if i0 <= ix <= i1
                    nray[ix, jyrv, kzrv] += 1
                    jray = nray[ix, jyrv, kzrv]
                    if jray > nray_wrk
                        error("Error in shift_rays!: nray > nray_wrk!")
                    end
                    copy_rays!(
                        rays,
                        (iray, ixrv, jyrv, kzrv),
                        (jray, ix, jyrv, kzrv),
                    )
                end
                rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
            end
        end
    end
end

function shift_rays!(state::State, direction::Y)
    (; sizezz, nzz, jo, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; ly, dy) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for kzrv in kz0:kz1, jyrv in (j0 - 1):(j1 + 1), ixrv in (i0 - 1):(i1 + 1)
        for iray in 1:nray[ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            jy = floor(Int, (yr - ly[1]) / dy) + j0 - jo

            if jy != jyrv
                if abs(jy - jyrv) > 1
                    error("Error in shift_rays!: abs(jy - jyrv) > 1!")
                end
                if j0 <= jy <= j1
                    nray[ixrv, jy, kzrv] += 1
                    jray = nray[ixrv, jy, kzrv]
                    if jray > nray_wrk
                        error("Error in shift_rays!: nray > nray_wrk!")
                    end
                    copy_rays!(
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

function shift_rays!(state::State, direction::Z)
    (; domain, grid) = state
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = domain
    (; nray_wrk, nray, rays) = state.wkb

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for kzrv in kz0:kz1, jyrv in (j0 - 1):(j1 + 1), ixrv in (i0 - 1):(i1 + 1)
        for iray in 1:nray[ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]
            kz = get_next_half_level(ixrv, jyrv, zr, domain, grid)

            if kz != kzrv
                if abs(kz - kzrv) > 1
                    error("Error in shift_rays!: abs(kz - kzrv) > 1!")
                end
                if k0 <= kz <= k1
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
                end
                rays.dens[iray, ixrv, jyrv, kzrv] = 0.0
            end
        end
    end

    return
end
