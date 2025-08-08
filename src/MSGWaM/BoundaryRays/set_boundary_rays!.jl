"""
```julia
set_boundary_rays!(state::State)
```

Enforce boundary conditions for ray volumes by dispatching to a test-case-specific method.

```julia
set_boundary_rays!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)
```

Enforce boundary conditions for ray volumes by dispatching to a WKB-mode-specific method.

```julia
set_boundary_rays!(state::State, wkb_mode::SteadyState)
```

Enforce horizontal boundary conditions for "ray volumes" in steady-state mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.sizex > 1` (`state.namelists.domain.sizey > 1`).

```julia
set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)
```

Enforce horizontal and vertical boundary conditions for ray volumes in single-column or multi-column mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.sizex > 1` (`state.namelists.domain.sizey > 1`).

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)
  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)
  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_boundary_rays!`](@ref)
"""
function set_boundary_rays! end

function set_boundary_rays!(state::State)
    (; testcase) = state.namelists.setting
    set_boundary_rays!(state, testcase)
    return
end

function set_boundary_rays!(state::State, testcase::AbstractTestCase)
    return
end

function set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_rays!(state, wkb_mode)
    return
end

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
