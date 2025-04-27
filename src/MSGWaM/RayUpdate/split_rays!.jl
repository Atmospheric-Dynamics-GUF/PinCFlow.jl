function split_rays!(state::State)
    (; testcase) = state.namelists.setting
    split_rays!(state, testcase)
    return
end

function split_rays!(state::State, testcase::AbstractTestCase)
    return
end

function split_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    split_rays!(state, wkb_mode)
    return
end

function split_rays!(state::State, wkb_mode::SteadyState)
    return
end

function split_rays!(state::State, wkb_mode::SingleColumn)
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray) = state.wkb

    @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        split_rays!(ix, jy, kz, state, Z())
    end

    @views nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_after = MPI.Allreduce(nray_after, +, comm)

    if master && nray_after > nray_before
        println("Number of ray volumes before splitting: ", nray_before)
        println("Number of ray volumes after splitting: ", nray_after)
        println("")
    end

    return
end

function split_rays!(state::State, wkb_mode::MultiColumn)
    (; sizex, sizey) = state.namelists.domain
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray) = state.wkb

    @views nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if sizex > 1
            split_rays!(ix, jy, kz, state, X())
        end

        if sizey > 1
            split_rays!(ix, jy, kz, state, Y())
        end

        split_rays!(ix, jy, kz, state, Z())
    end

    @views nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_after = MPI.Allreduce(nray_after, +, comm)

    if master && nray_after > nray_before
        println("Number of ray volumes before splitting: ", nray_before)
        println("Number of ray volumes after splitting: ", nray_after)
        println("")
    end

    return
end

function split_rays!(
    ix::Integer,
    jy::Integer,
    kz::Integer,
    state::State,
    axis::X,
)
    (; dx) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    nrlc = nray[ix, jy, kz]
    for iray in 1:nray[ix, jy, kz]
        xr = rays.x[iray, ix, jy, kz]
        dxr = rays.dxray[iray, ix, jy, kz]

        if dxr > dx
            nrlc += 1

            rays.dxray[iray, ix, jy, kz] = 0.5 * dxr

            copy_rays!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))

            rays.x[iray, ix, jy, kz] = xr - 0.25 * dxr
            rays.x[nrlc, ix, jy, kz] = xr + 0.25 * dxr
        end
    end

    if nrlc > nray[ix, jy, kz]
        nray[ix, jy, kz] = nrlc

        if nray[ix, jy, kz] > nray_wrk
            error(
                "Error in split_rays!: nray",
                [ix, jy, iz],
                " > nray_wrk = ",
                nray_wrk,
            )
        end
    end

    return
end

function split_rays!(
    ix::Integer,
    jy::Integer,
    kz::Integer,
    state::State,
    axis::Y,
)
    (; dy) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    nrlc = nray[ix, jy, kz]
    for iray in 1:nray[ix, jy, kz]
        yr = rays.y[iray, ix, jy, kz]
        dyr = rays.dyray[iray, ix, jy, kz]

        if dyr > dy
            nrlc += 1

            rays.dyray[iray, ix, jy, kz] = 0.5 * dyr

            copy_rays!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))

            rays.y[iray, ix, jy, kz] = yr - 0.25 * dyr
            rays.y[nrlc, ix, jy, kz] = yr + 0.25 * dyr
        end
    end

    if nrlc > nray[ix, jy, kz]
        nray[ix, jy, kz] = nrlc

        if nray[ix, jy, kz] > nray_wrk
            error(
                "Error in split_rays!: nray",
                [ix, jy, iz],
                " > nray_wrk = ",
                nray_wrk,
            )
        end
    end

    return
end

function split_rays!(
    ix::Integer,
    jy::Integer,
    kz::Integer,
    state::State,
    axis::Z,
)
    (; domain, grid) = state
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, dz, jac) = grid
    (; nray_wrk, nray, rays) = state.wkb

    nrlc = nray[ix, jy, kz]
    for iray in 1:nray[ix, jy, kz]
        xr = rays.x[iray, ix, jy, kz]
        yr = rays.y[iray, ix, jy, kz]
        zr = rays.z[iray, ix, jy, kz]

        dzr = rays.dzray[iray, ix, jy, kz]

        ixrv = floor(Int, (xr - lx[1]) / dx) + i0 - io
        jyrv = floor(Int, (yr - ly[1]) / dy) + j0 - jo
        kzrvd = get_next_half_level(ixrv, jyrv, zr - 0.5 * dzr, domain, grid)
        kzrvu = get_next_half_level(ixrv, jyrv, zr + 0.5 * dzr, domain, grid)

        dzmin = dz
        for kzrv in kzrvd:kzrvu
            dzmin = min(dzmin, jac[ixrv, jyrv, kzrv] * dz)
        end

        if dzr > dzmin
            factor = ceil(Int, dzr / dzmin)
            rays.z[iray, ix, jy, kz] = zr + 0.5 * (1 / factor - 1) * dzr
            rays.dzray[iray, ix, jy, kz] = dzr / factor
            for jray in (nrlc + 1):(nrlc + factor - 1)
                copy_rays!(rays, (iray, ix, jy, kz), (jray, ix, jy, kz))
                rays.z[jray, ix, jy, kz] += (jray - nrlc) * dzr / factor
            end
            nrlc += factor - 1
        end
    end

    if nrlc > nray[ix, jy, kz]
        nray[ix, jy, kz] = nrlc

        if nray[ix, jy, kz] > nray_wrk
            error(
                "Error in split_rays!: nray",
                [ix, jy, iz],
                " > nray_wrk = ",
                nray_wrk,
            )
        end
    end

    return
end
