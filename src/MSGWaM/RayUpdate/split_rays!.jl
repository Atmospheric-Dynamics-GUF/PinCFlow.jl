"""
```julia
split_rays!(state::State)
```

Entry point for ray splitting operations. Delegates to the appropriate method based on the test case specified in the state's settings.

# Arguments

  - `state::State`: The state object containing model data and configuration
"""
function split_rays!(state::State)
    (; testcase) = state.namelists.setting
    split_rays!(state, testcase)
    return
end

"""
```julia
split_rays!(state::State, testcase::AbstractTestCase)
```

Default implementation for non-WKB test cases. Does nothing as ray splitting is not required for standard test cases.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `testcase::AbstractTestCase`: The test case specification
"""
function split_rays!(state::State, testcase::AbstractTestCase)
    return
end

"""
```julia
split_rays!(state::State, testcase::AbstractWKBTestCase)
```

Implementation for WKB test cases. Delegates to the appropriate method based on
the WKB mode specified in the state's configuration.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `testcase::AbstractWKBTestCase`: The WKB test case specification
"""
function split_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    split_rays!(state, wkb_mode)
    return
end

"""
```julia
split_rays!(state::State, wkb_mode::SteadyState)
```

Implementation for steady state WKB mode. Does nothing as rays remain stationary
and do not require splitting.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `wkb_mode::SteadyState`: The steady state WKB mode specification
"""
function split_rays!(state::State, wkb_mode::SteadyState)
    return
end

"""
```julia
split_rays!(state::State, wkb_mode::SingleColumn)
```

Implementation for single column WKB mode. Splits rays in the vertical (Z) direction only
when their extent becomes too large relative to the grid spacing.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `wkb_mode::SingleColumn`: The single column WKB mode specification

# Implementation

  - Counts total ray volumes before and after splitting across all MPI processes
  - Performs vertical ray splitting for each grid cell
  - Reports splitting statistics if rays were actually split
  - Only the master process prints output to avoid duplicate messages
"""
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

"""
```julia
split_rays!(state::State, wkb_mode::MultiColumn)
```

Implementation for multi-column WKB mode. Splits rays in X, Y, and Z directions
as appropriate based on domain dimensions.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `wkb_mode::MultiColumn`: The multi-column WKB mode specification

# Implementation

  - Counts total ray volumes before and after splitting across all MPI processes

  - For each grid cell, performs splitting in active dimensions:

      + X-direction splitting if domain size in X > 1
      + Y-direction splitting if domain size in Y > 1
      + Z-direction splitting (always performed)
  - Reports splitting statistics if rays were actually split
  - Only the master process prints output to avoid duplicate messages
"""
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

"""
```julia
split_rays!(ix::Integer, jy::Integer, kz::Integer, state::State, axis::X)
```

Split rays in the X (zonal) direction for a specific grid cell when ray extent
exceeds the grid spacing.

# Arguments

  - `ix::Integer`: Grid cell index in X direction
  - `jy::Integer`: Grid cell index in Y direction
  - `kz::Integer`: Grid cell index in Z direction
  - `state::State`: The state object containing model data and configuration
  - `axis::X`: Dispatch type indicating X-direction splitting

# Implementation

  - Examines each ray's X-direction extent (`dxray`)

  - If `dxray > dx` (grid spacing), splits the ray into two parts:

      + Both parts have half the original extent (`0.5 * dxray`)
      + Original ray position shifted to `xr - 0.25 * dxray`
      + New ray position set to `xr + 0.25 * dxray`
  - Updates ray count and performs bounds checking
  - Copies all ray properties to the new ray volume

# Error Handling

  - Throws error if total ray count exceeds `nray_wrk` limit
"""
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

"""
```julia
split_rays!(ix::Integer, jy::Integer, kz::Integer, state::State, axis::Y)
```

Split rays in the Y (meridional) direction for a specific grid cell when ray extent
exceeds the grid spacing.

# Arguments

  - `ix::Integer`: Grid cell index in X direction
  - `jy::Integer`: Grid cell index in Y direction
  - `kz::Integer`: Grid cell index in Z direction
  - `state::State`: The state object containing model data and configuration
  - `axis::Y`: Dispatch type indicating Y-direction splitting

# Implementation

  - Examines each ray's Y-direction extent (`dyray`)

  - If `dyray > dy` (grid spacing), splits the ray into two parts:

      + Both parts have half the original extent (`0.5 * dyray`)
      + Original ray position shifted to `yr - 0.25 * dyray`
      + New ray position set to `yr + 0.25 * dyray`
  - Updates ray count and performs bounds checking
  - Copies all ray properties to the new ray volume

# Error Handling

  - Throws error if total ray count exceeds `nray_wrk` limit
"""
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

"""
```julia
split_rays!(ix::Integer, jy::Integer, kz::Integer, state::State, axis::Z)
```

Split rays in the Z (vertical) direction for a specific grid cell when ray extent
exceeds the minimum local grid spacing.

# Arguments

  - `ix::Integer`: Grid cell index in X direction
  - `jy::Integer`: Grid cell index in Y direction
  - `kz::Integer`: Grid cell index in Z direction
  - `state::State`: The state object containing model data and configuration
  - `axis::Z`: Dispatch type indicating Z-direction splitting

# Implementation

  - Determines vertical extent of ray (`dzray`) and its spatial range

  - Finds minimum grid spacing (`dzmin`) over all levels the ray spans
  - Uses `get_next_half_level` to determine vertical level boundaries
  - Accounts for grid stretching via Jacobian factor
  - If `dzray > dzmin`, splits ray into multiple parts:

      + Number of parts: `factor = ceil(dzray / dzmin)`
      + Each part has extent `dzray / factor`
      + Parts are distributed vertically with spacing `dzray / factor`
  - Updates ray count and performs bounds checking

# Error Handling

  - Throws error if total ray count exceeds `nray_wrk` limit

# Notes

  - More complex than X/Y splitting due to vertical grid stretching
  - Can split into more than 2 parts if necessary
  - Respects terrain-following coordinate system
"""
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
            rays.z[iray, ix, jy, kz] = zr + (1 / factor - 1) * dzr / 2
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
