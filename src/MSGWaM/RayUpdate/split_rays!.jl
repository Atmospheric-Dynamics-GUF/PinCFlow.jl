"""
```julia
split_rays!(state::State)
```

Split ray volumes that have become larger than the local grid cell by dispatching to a test-case-specific method.

```julia
split_rays!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
split_rays!(state::State, testcase::AbstractWKBTestCase)
```

Split ray volumes that have become larger than the local grid cell by dispatching to a WKB-mode-specific method.

```julia
split_rays!(state::State, wkb_mode::SteadyState)
```

Return for steady-state mode.

```julia
split_rays!(state::State, wkb_mode::SingleColumn)
```

Split ray volumes which have a vertical extent larger than the local vertical grid spacing.

```julia
split_rays!(state::State, wkb_mode::MultiColumn)
```

In each dimension of physical space, split ray volumes which have an extent larger than the local grid spacing.

The splitting is performed sequentially, such that a ray volume with extents that are all between once and twice as large as allowed is split into exactly eight smaller ray volumes (all of which have the same size).

```julia
split_rays!(ix::Integer, jy::Integer, kz::Integer, state::State, axis::X)
```

In the grid cell specified by `ix`, `jy` and `kz`, split ray volumes with ``\\Delta x_\\alpha > \\Delta \\widehat{x}``.

The number of splits is the result of ceiling division of ``\\Delta x_\\alpha`` by ``\\Delta \\widehat{x}``. Each split is carried out by adjusting the position and extent of the ray volume, copying it and changing the position of the copy appropriately.

```julia
split_rays!(ix::Integer, jy::Integer, kz::Integer, state::State, axis::Y)
```

In the grid cell specified by `ix`, `jy` and `kz`, split ray volumes with ``\\Delta y_\\alpha > \\Delta \\widehat{y}``.

The splitting is analogous to that in ``\\widehat{x}``.

```julia
split_rays!(ix::Integer, jy::Integer, kz::Integer, state::State, axis::Z)
```

In the grid cell specified by `ix`, `jy` and `kz`, split ray volumes with ``\\Delta z_\\alpha > J_{\\min} \\Delta \\widehat{z}``, with ``J_{\\min}`` being the minimum value of the Jacobian in all grid cells that are at least partially covered by the ray volume (at its true horizontal position on the grid).

The splitting is analogous to that in ``\\widehat{x}`` and ``\\widehat{y}``.

# Arguments

  - `state`: Model state.

  - `testcase`: Test case on which the current simulation is based.

  - `wkb_mode`: Approximations used by MSGWaM.

  - `ix`: Grid-cell index in ``\\widehat{x}``-direction

  - `jy`: Grid-cell index in ``\\widehat{y}``-direction

  - `kz`: Grid-cell index in ``\\widehat{z}``-direction

  - `axis`: Axis perpendicular to the split.

# See also

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function split_rays! end

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

    @ivy nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    @ivy for kz in k0:k1, jy in j0:j1, ix in i0:i1
        split_rays!(ix, jy, kz, state, Z())
    end

    @ivy nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
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

    @ivy nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    @ivy for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if sizex > 1
            split_rays!(ix, jy, kz, state, X())
        end

        if sizey > 1
            split_rays!(ix, jy, kz, state, Y())
        end

        split_rays!(ix, jy, kz, state, Z())
    end

    @ivy nray_after = sum(nray[i0:i1, j0:j1, k0:k1])
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

    @ivy nrlc = nray[ix, jy, kz]
    @ivy for iray in 1:nray[ix, jy, kz]
        xr = rays.x[iray, ix, jy, kz]
        dxr = rays.dxray[iray, ix, jy, kz]

        if dxr > dx
            factor = ceil(Int, dxr / dx)
            rays.x[iray, ix, jy, kz] = xr + (1 / factor - 1) * dxr / 2
            rays.dxray[iray, ix, jy, kz] = dxr / factor
            for jray in (nrlc + 1):(nrlc + factor - 1)
                copy_rays!(rays, (iray, ix, jy, kz), (jray, ix, jy, kz))
                rays.x[jray, ix, jy, kz] += (jray - nrlc) * dxr / factor
            end
            nrlc += factor - 1
        end
    end

    @ivy if nrlc > nray[ix, jy, kz]
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

    @ivy nrlc = nray[ix, jy, kz]
    @ivy for iray in 1:nray[ix, jy, kz]
        yr = rays.y[iray, ix, jy, kz]
        dyr = rays.dyray[iray, ix, jy, kz]

        if dyr > dy
            factor = ceil(Int, dyr / dy)
            rays.y[iray, ix, jy, kz] = yr + (1 / factor - 1) * dyr / 2
            rays.dyray[iray, ix, jy, kz] = dyr / factor
            for jray in (nrlc + 1):(nrlc + factor - 1)
                copy_rays!(rays, (iray, ix, jy, kz), (jray, ix, jy, kz))
                rays.y[jray, ix, jy, kz] += (jray - nrlc) * dyr / factor
            end
            nrlc += factor - 1
        end
    end

    @ivy if nrlc > nray[ix, jy, kz]
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

    @ivy nrlc = nray[ix, jy, kz]
    @ivy for iray in 1:nray[ix, jy, kz]
        xr = rays.x[iray, ix, jy, kz]
        yr = rays.y[iray, ix, jy, kz]
        zr = rays.z[iray, ix, jy, kz]

        dzr = rays.dzray[iray, ix, jy, kz]

        ixrv = floor(Int, (xr + lx / 2) / dx) + i0 - io
        jyrv = floor(Int, (yr + ly / 2) / dy) + j0 - jo
        kzrvd = get_next_half_level(ixrv, jyrv, zr - 0.5 * dzr, domain, grid)
        kzrvu = get_next_half_level(ixrv, jyrv, zr + 0.5 * dzr, domain, grid)

        dzmin = dz
        for kzrv in kzrvd:kzrvu
            dzmin = min(dzmin, jac[ixrv, jyrv, kzrv] * dz)
        end

        if dzr > dzmin
            factor = ceil(Int, dzr / dzmin)
            rays.z[iray, ix, jy, kz] = zr + (1 / factor - 1) * dzr / 2
            rays.dzray[iray, ix, jy, kz] = dzr / factor
            for jray in (nrlc + 1):(nrlc + factor - 1)
                copy_rays!(rays, (iray, ix, jy, kz), (jray, ix, jy, kz))
                rays.z[jray, ix, jy, kz] += (jray - nrlc) * dzr / factor
            end
            nrlc += factor - 1
        end
    end

    @ivy if nrlc > nray[ix, jy, kz]
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
