"""
```julia
split_rays!(state::State)
```

Split ray volumes that have become larger than the local grid cell by dispatching to a WKB-mode-specific method.

```julia
split_rays!(state::State, wkb_mode::Union{NoWKB, SteadyState})
```

Return for configurations without WKB / with steady-state WKB.

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
split_rays!(i::Integer, j::Integer, k::Integer, state::State, axis::X)
```

In the grid cell specified by ``\\left(i, j, k\\right)``, split ray volumes with ``\\Delta x_r > \\Delta \\widehat{x}``.

The number of splits is the result of ceiling division of ``\\Delta x_r`` by ``\\Delta \\widehat{x}``. Each split is carried out by adjusting the position and extent of the ray volume, copying it and changing the position of the copy appropriately.

```julia
split_rays!(i::Integer, j::Integer, k::Integer, state::State, axis::Y)
```

In the grid cell specified by ``\\left(i, j, k\\right)``, split ray volumes with ``\\Delta y_r > \\Delta \\widehat{y}``.

The splitting is analogous to that in ``\\widehat{x}``.

```julia
split_rays!(i::Integer, j::Integer, k::Integer, state::State, axis::Z)
```

In the grid cell specified by ``\\left(i, j, k\\right)``, split ray volumes with ``\\Delta z_r > J_{\\min} \\Delta \\widehat{z}``, with ``J_{\\min}`` being the minimum value of the Jacobian in all grid cells that are at least partially covered by the ray volume (at its true horizontal position on the grid).

The splitting is analogous to that in ``\\widehat{x}`` and ``\\widehat{y}``.

# Arguments

  - `state`: Model state.

  - `wkb_mode`: Approximations used by MS-GWaM.

  - `i`: Grid-cell index in ``\\widehat{x}``-direction

  - `j`: Grid-cell index in ``\\widehat{y}``-direction

  - `k`: Grid-cell index in ``\\widehat{z}``-direction

  - `axis`: Axis perpendicular to the split.

# See also

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
"""
function split_rays! end

function split_rays!(state::State)
    (; wkb_mode) = state.namelists.wkb
    split_rays!(state, wkb_mode)
    return
end

function split_rays!(state::State, wkb_mode::Union{NoWKB, SteadyState})
    return
end

function split_rays!(state::State, wkb_mode::SingleColumn)
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray) = state.wkb

    @ivy nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        split_rays!(i, j, k, state, Z())
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
    (; x_size, y_size) = state.namelists.domain
    (; comm, master, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray) = state.wkb

    @ivy nray_before = sum(nray[i0:i1, j0:j1, k0:k1])
    nray_before = MPI.Allreduce(nray_before, +, comm)

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        if x_size > 1
            split_rays!(i, j, k, state, X())
        end

        if y_size > 1
            split_rays!(i, j, k, state, Y())
        end

        split_rays!(i, j, k, state, Z())
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

function split_rays!(i::Integer, j::Integer, k::Integer, state::State, axis::X)
    (; dx) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    @ivy local_count = nray[i, j, k]
    @ivy for r in 1:nray[i, j, k]
        xr = rays.x[r, i, j, k]
        dxr = rays.dxray[r, i, j, k]

        if dxr > dx
            factor = ceil(Int, dxr / dx)
            rays.x[r, i, j, k] = xr + (1 / factor - 1) * dxr / 2
            rays.dxray[r, i, j, k] = dxr / factor
            for rray in (local_count + 1):(local_count + factor - 1)
                copy_rays!(rays, r => rray, i => i, j => j, k => k)
                rays.x[rray, i, j, k] += (rray - local_count) * dxr / factor
            end
            local_count += factor - 1
        end
    end

    @ivy if local_count > nray[i, j, k]
        nray[i, j, k] = local_count

        if nray[i, j, k] > nray_wrk
            error(
                "Error in split_rays!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
        end
    end

    return
end

function split_rays!(i::Integer, j::Integer, k::Integer, state::State, axis::Y)
    (; dy) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    @ivy local_count = nray[i, j, k]
    @ivy for r in 1:nray[i, j, k]
        yr = rays.y[r, i, j, k]
        dyr = rays.dyray[r, i, j, k]

        if dyr > dy
            factor = ceil(Int, dyr / dy)
            rays.y[r, i, j, k] = yr + (1 / factor - 1) * dyr / 2
            rays.dyray[r, i, j, k] = dyr / factor
            for rray in (local_count + 1):(local_count + factor - 1)
                copy_rays!(rays, r => rray, i => i, j => j, k => k)
                rays.y[rray, i, j, k] += (rray - local_count) * dyr / factor
            end
            local_count += factor - 1
        end
    end

    @ivy if local_count > nray[i, j, k]
        nray[i, j, k] = local_count

        if nray[i, j, k] > nray_wrk
            error(
                "Error in split_rays!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
        end
    end

    return
end

function split_rays!(i::Integer, j::Integer, k::Integer, state::State, axis::Z)
    (; domain, grid) = state
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, dz, jac) = grid
    (; nray_wrk, nray, rays) = state.wkb

    @ivy local_count = nray[i, j, k]
    @ivy for r in 1:nray[i, j, k]
        xr = rays.x[r, i, j, k]
        yr = rays.y[r, i, j, k]
        zr = rays.z[r, i, j, k]

        dzr = rays.dzray[r, i, j, k]

        iray = floor(Int, (xr + lx / 2) / dx) + i0 - io
        jray = floor(Int, (yr + ly / 2) / dy) + j0 - jo
        kmin = get_next_half_level(iray, jray, zr - 0.5 * dzr, state)
        kmax = get_next_half_level(iray, jray, zr + 0.5 * dzr, state)

        dzmin = dz
        for kray in kmin:kmax
            dzmin = min(dzmin, jac[iray, jray, kray] * dz)
        end

        if dzr > dzmin
            factor = ceil(Int, dzr / dzmin)
            rays.z[r, i, j, k] = zr + (1 / factor - 1) * dzr / 2
            rays.dzray[r, i, j, k] = dzr / factor
            for rray in (local_count + 1):(local_count + factor - 1)
                copy_rays!(rays, r => rray, i => i, j => j, k => k)
                rays.z[rray, i, j, k] += (rray - local_count) * dzr / factor
            end
            local_count += factor - 1
        end
    end

    @ivy if local_count > nray[i, j, k]
        nray[i, j, k] = local_count

        if nray[i, j, k] > nray_wrk
            error(
                "Error in split_rays!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
        end
    end

    return
end
