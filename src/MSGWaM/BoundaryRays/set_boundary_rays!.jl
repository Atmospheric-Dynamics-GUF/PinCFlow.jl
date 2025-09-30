"""
```julia
set_boundary_rays!(state::State)
```

Enforce boundary conditions for ray volumes by dispatching to a test-case-specific method.

```julia
set_boundary_rays!(state::State, test_case::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
set_boundary_rays!(state::State, test_case::AbstractWKBTestCase)
```

Enforce boundary conditions for ray volumes by dispatching to a WKB-mode-specific method.

```julia
set_boundary_rays!(state::State, wkb_mode::SteadyState)
```

Enforce horizontal boundary conditions for "ray volumes" in steady-state mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.ndx > 1` (`state.namelists.domain.ndy > 1`).

```julia
set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)
```

Enforce horizontal and vertical boundary conditions for ray volumes in single-column or multi-column mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.ndx > 1` (`state.namelists.domain.ndy > 1`).

# Arguments

  - `state`: Model state.

  - `test_case`: Test case on which the current simulation is based.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_boundary_rays!`](@ref)
"""
function set_boundary_rays! end

function set_boundary_rays!(state::State)
    (; test_case) = state.namelists.setting
    set_boundary_rays!(state, test_case)
    return
end

function set_boundary_rays!(state::State, test_case::AbstractTestCase)
    return
end

function set_boundary_rays!(state::State, test_case::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_rays!(state, wkb_mode)
    return
end

function set_boundary_rays!(state::State, wkb_mode::SteadyState)
    (; ndx, ndy) = state.namelists.domain

    if ndx > 1
        set_zonal_boundary_rays!(state)
    end
    if ndy > 1
        set_meridional_boundary_rays!(state)
    end

    return
end

function set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)
    (; ndx, ndy) = state.namelists.domain

    if ndx > 1
        set_zonal_boundary_rays!(state)
    end
    if ndy > 1
        set_meridional_boundary_rays!(state)
    end
    set_vertical_boundary_rays!(state)

    return
end
