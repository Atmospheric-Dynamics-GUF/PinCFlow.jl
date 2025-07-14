"""
```julia
set_boundary_rays!(state::State)
```

Enforce boundary conditions for ray volumes based on the test case.

Dispatches to the appropriate methods depending on whether or not the simulation uses MS-GWaM (i.e. whether or not a WKB test case is simulated).

# Arguments

  - `state`: Model state.
"""
function set_boundary_rays!(state::State)
    (; testcase) = state.namelists.setting
    set_boundary_rays!(state, testcase)
    return
end

"""
```julia
set_boundary_rays!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.
"""
function set_boundary_rays!(state::State, testcase::AbstractTestCase)
    return
end

"""
```julia
set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)
```

Enforce boundary conditions for ray volumes based on the WKB mode.

Dispatches to specific methods depending on the WKB mode.

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.
"""
function set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_rays!(state, wkb_mode)
    return
end

"""
```julia
set_boundary_rays!(state::State, wkb_mode::SteadyState)
```

Enforce horizontal boundary conditions for "ray volumes" in steady-state mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.sizex > 1` (`state.namelists.domain.sizey > 1`).

# Arguments

  - `state`: Model state.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)
  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)
"""
function set_boundary_rays!(state::State, wkb_mode::SteadyState)
    (; sizex, sizey) = state.namelists.domain

    if sizex > 1
        set_zonal_boundary_rays!(state)
    end
    if sizey > 1
        set_meridional_boundary_rays!(state)
    end

    return
end

"""
```julia
set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)
```

Enforce horizontal and vertical boundary conditions for ray volumes in single-column or multi-column mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.sizex > 1` (`state.namelists.domain.sizey > 1`).

# Arguments

  - `state`: Model state.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)
  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)
  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_boundary_rays!`](@ref)
"""
function set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)
    (; sizex, sizey) = state.namelists.domain
    (; zboundaries) = state.namelists.setting

    if sizex > 1
        set_zonal_boundary_rays!(state)
    end
    if sizey > 1
        set_meridional_boundary_rays!(state)
    end
    set_vertical_boundary_rays!(state, zboundaries)

    return
end
