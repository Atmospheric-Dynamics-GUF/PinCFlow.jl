"""
```julia
shift_rays!(state::State)
```

Shift the array positions of ray volumes such that they are attributed to the correct grid cells by dispatching to a test-case-specific method.

```julia
shift_rays!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
shift_rays!(state::State, testcase::AbstractWKBTestCase)
```

Shift the array positions of ray volumes such that they are attributed to the correct grid cells by dispatching to a WKB-mode-specific method.

```julia
shift_rays!(state::State, wkb_mode::SteadyState)
```

Return for steady-state mode.

```julia
shift_rays!(state::State, wkb_mode::SingleColumn)
```

Shift the vertical array positions of ray volumes such that they are attributed to the correct grid cells.

This method enforces the vertical boundary conditions (via `set_vertical_boundary_rays!`), checks if ray volumes need to be shifted and, if they do, copies them to the correct grid cells and marks them for removal (by dispatching to the appropriate method). A second call of `set_vertical_boundary_rays!` ensures that ray volumes that have moved across MPI processes are included in the appropriate halo cells. Finally, the gaps that were created by marking ray volumes for removal are filled (via `remove_rays!`).

```julia
shift_rays!(state::State, wkb_mode::MultiColumn)
```

Shift the array positions of ray volumes such that they are attributed to the correct grid cells.

For each dimension in physical space (with more than one grid point), this method performs the corresponding equivalent of the algorithm that is implemented in the method for single-column mode.

```julia
shift_rays!(state::State, direction::X)
```

For each ray volume, check if it is attributed to the correct position in ``\\widehat{x}`` and, if it is not, create a copy that is and mark the original for removal.

Ray volumes that should be attributed to a halo cell are marked for removal but not copied, since the copies are created from the corresponding halo cell in the adjacent MPI process.

```julia
shift_rays!(state::State, direction::Y)
```

For each ray volume, check if it is attributed to the correct position in ``\\widehat{y}`` and, if it is not, create a copy that is and mark the original for removal.

Ray volumes in halo cells are treated in the same way as in the method for shifting in ``\\widehat{x}``.

```julia
shift_rays!(state::State, direction::Z)
```

For each ray volume, check if it is attributed to the correct position in ``\\widehat{z}`` and, if it is not, create a copy that is and mark the original for removal.

Ray volumes in halo cells are treated in the same way as in the methods for shifting in ``\\widehat{x}`` and ``\\widehat{z}``.

# Arguments

  - `state`: Model state.

  - `testcase`: Test case on which the current simulation is based.

  - `wkb_mode`: Approximations used by MSGWaM.

  - `direction`: Shift direction.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.remove_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.check_rays`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function shift_rays! end

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
    set_vertical_boundary_rays!(state)
    shift_rays!(state, Z())
    set_vertical_boundary_rays!(state)
    remove_rays!(state)

    check_rays(state)

    return
end

function shift_rays!(state::State, wkb_mode::MultiColumn)
    (; sizex, sizey) = state.namelists.domain

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

    set_vertical_boundary_rays!(state)
    shift_rays!(state, Z())
    set_vertical_boundary_rays!(state)
    remove_rays!(state)

    check_rays(state)

    return
end

function shift_rays!(state::State, direction::X)
    (; sizezz, nzz, io, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, dx) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    ks = ko == 0 ? k0 : k0 - 1
    ke = ko + nzz == sizezz ? k1 : k1 + 1

    @ivy for k in ks:ke, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
        for r in 1:nray[i, j, k]
            xr = rays.x[r, i, j, k]
            iray = floor(Int, (xr + lx / 2) / dx) + i0 - io

            if iray != i
                if abs(iray - i) > 1
                    error("Error in shift_rays!: abs(iray - i) > 1!")
                end
                if i0 <= iray <= i1
                    nray[iray, j, k] += 1
                    rray = nray[iray, j, k]
                    if rray > nray_wrk
                        error("Error in shift_rays!: nray > nray_wrk!")
                    end
                    copy_rays!(rays, r => rray, i => iray, j => j, k => k)
                end
                rays.dens[r, i, j, k] = 0.0
            end
        end
    end

    return
end

function shift_rays!(state::State, direction::Y)
    (; sizezz, nzz, jo, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; ly, dy) = state.grid
    (; nray_wrk, nray, rays) = state.wkb

    ks = ko == 0 ? k0 : k0 - 1
    ke = ko + nzz == sizezz ? k1 : k1 + 1

    @ivy for k in ks:ke, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
        for r in 1:nray[i, j, k]
            yr = rays.y[r, i, j, k]
            jray = floor(Int, (yr + ly / 2) / dy) + j0 - jo

            if jray != j
                if abs(jray - j) > 1
                    error("Error in shift_rays!: abs(jray - j) > 1!")
                end
                if j0 <= jray <= j1
                    nray[i, jray, k] += 1
                    rray = nray[i, jray, k]
                    if rray > nray_wrk
                        error("Error in shift_rays!: nray > nray_wrk!")
                    end
                    copy_rays!(rays, r => rray, i => i, j => jray, k => k)
                end
                rays.dens[r, i, j, k] = 0.0
            end
        end
    end

    return
end

function shift_rays!(state::State, direction::Z)
    (; domain, grid) = state
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = domain
    (; nray_wrk, nray, rays) = state.wkb

    ks = ko == 0 ? k0 : k0 - 1
    ke = ko + nzz == sizezz ? k1 : k1 + 1

    @ivy for k in ks:ke, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
        for r in 1:nray[i, j, k]
            zr = rays.z[r, i, j, k]
            kray = get_next_half_level(i, j, zr, state)

            if kray != k
                if abs(kray - k) > 1
                    error("Error in shift_rays!: abs(kray - k) > 1!")
                end
                if k0 <= kray <= k1
                    nray[i, j, kray] += 1
                    rray = nray[i, j, kray]
                    if rray > nray_wrk
                        error("Error in shift_rays!: nray > nray_wrk!")
                    end
                    copy_rays!(rays, r => rray, i => i, j => j, k => kray)
                end
                rays.dens[r, i, j, k] = 0.0
            end
        end
    end

    return
end
