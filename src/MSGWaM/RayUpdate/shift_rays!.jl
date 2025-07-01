"""
```julia
shift_rays!(state::State)
```

Entry point for ray shifting operations. Delegates to the appropriate method based on
the test case specified in the state's settings.

# Arguments

  - `state::State`: The state object containing model data and configuration
"""
function shift_rays!(state::State)
    (; testcase) = state.namelists.setting
    shift_rays!(state, testcase)
    return
end

"""
```julia
shift_rays!(state::State, testcase::AbstractTestCase)
```

Default implementation for non-WKB test cases. Does nothing.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `testcase::AbstractTestCase`: The test case specification
"""
function shift_rays!(state::State, testcase::AbstractTestCase)
    return
end

"""
```julia
shift_rays!(state::State, testcase::AbstractWKBTestCase)
```

Implementation for WKB test cases. Delegates to the appropriate method based on
the WKB mode specified in the state's configuration.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `testcase::AbstractWKBTestCase`: The WKB test case specification
"""
function shift_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    shift_rays!(state, wkb_mode)
    return
end

"""
```julia
shift_rays!(state::State, wkb_mode::SteadyState)
```

Implementation for steady state WKB mode. Does nothing as rays remain stationary.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `wkb_mode::SteadyState`: The steady state WKB mode specification
"""
function shift_rays!(state::State, wkb_mode::SteadyState)
    return
end

"""
```julia
shift_rays!(state::State, wkb_mode::SingleColumn)
```

Implementation for single column WKB mode. Shifts rays in the vertical (Z) direction only.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `wkb_mode::SingleColumn`: The single column WKB mode specification

# Workflow

 1. Sets vertical boundary rays
 2. Shifts rays in Z direction
 3. Re-establishes vertical boundary rays
 4. Removes invalid rays
 5. Performs ray validation checks
"""
function shift_rays!(state::State, wkb_mode::SingleColumn)
    (; zboundaries) = state.namelists.setting

    set_vertical_boundary_rays!(state, zboundaries)
    shift_rays!(state, Z())
    set_vertical_boundary_rays!(state, zboundaries)
    remove_rays!(state)

    check_rays(state)

    return
end

"""
```julia
shift_rays!(state::State, wkb_mode::MultiColumn)
```

Implementation for multi-column WKB mode. Shifts rays in X, Y, and Z directions
as appropriate based on domain dimensions.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `wkb_mode::MultiColumn`: The multi-column WKB mode specification

# Workflow

 1. If X-dimension > 1: Sets zonal boundary rays, shifts rays in X direction,
    re-establishes zonal boundary rays, and removes invalid rays
 2. If Y-dimension > 1: Sets meridional boundary rays, shifts rays in Y direction,
    re-establishes meridional boundary rays, and removes invalid rays
 3. Sets vertical boundary rays, shifts rays in Z direction, re-establishes vertical
    boundary rays, and removes invalid rays
 4. Performs ray validation checks
"""
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

"""
```julia
shift_rays!(state::State, direction::X)
```

Shifts rays in the X (zonal) direction when their positions no longer correspond
to their assigned grid cell.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `direction::X`: Dispatch type indicating X-direction shifting

# Implementation

  - Examines each ray's X position
  - If a ray has moved to a new grid cell, transfers it to the appropriate cell
  - Sets density to zero in the original location to mark the original ray as invalid
  - Error checks to ensure rays don't move more than one grid cell at a time
"""
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

    return
end

"""
```julia
shift_rays!(state::State, direction::Y)
```

Shifts rays in the Y (meridional) direction when their positions no longer correspond
to their assigned grid cell.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `direction::Y`: Dispatch type indicating Y-direction shifting

# Implementation

  - Examines each ray's Y position
  - If a ray has moved to a new grid cell, transfers it to the appropriate cell
  - Sets density to zero in the original location to mark the original ray as invalid
  - Error checks to ensure rays don't move more than one grid cell at a time
"""
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

    return
end

"""
```julia
shift_rays!(state::State, direction::Z)
```

Shifts rays in the Z (vertical) direction when their positions no longer correspond
to their assigned grid cell.

# Arguments

  - `state::State`: The state object containing model data and configuration
  - `direction::Z`: Dispatch type indicating Z-direction shifting

# Implementation

  - Examines each ray's Z position
  - Uses `get_next_half_level` to determine the appropriate vertical level
  - If a ray has moved to a new grid cell, transfers it to the appropriate cell
  - Sets density to zero in the original location to mark the original ray as invalid
  - Enforces domain boundaries for vertical movement
"""
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
