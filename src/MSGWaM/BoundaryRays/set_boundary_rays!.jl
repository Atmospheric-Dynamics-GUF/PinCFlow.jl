"""
```julia
set_boundary_rays!(state::State)
```

Enforce boundary conditions for ray volumes by dispatching to a WKB-mode-specific method.

```julia
set_boundary_rays!(state::State, wkb_mode::NoWKB)
```

Return for non-WKB configurations.

```julia
set_boundary_rays!(state::State, wkb_mode::SteadyState)
```

Enforce horizontal boundary conditions for "ray volumes" in steady-state mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.x_size > 1` (`state.namelists.domain.y_size > 1`).

```julia
set_boundary_rays!(state::State, wkb_mode::Union{SingleColumn, MultiColumn})
```

Enforce horizontal and vertical boundary conditions for ray volumes in single-column or multi-column mode.

Zonal (meridional) boundary conditions are only enforced if `state.namelists.domain.x_size > 1` (`state.namelists.domain.y_size > 1`).

# Arguments

  - `state`: Model state.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_boundary_rays!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_boundary_rays!`](@ref)
"""
function set_boundary_rays! end

function set_boundary_rays!(state::State)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_rays!(state, wkb_mode)
    return
end

function set_boundary_rays!(state::State, wkb_mode::NoWKB)
    return
end

function set_boundary_rays!(state::State, wkb_mode::SteadyState)
    (; x_size, y_size) = state.namelists.domain

    if x_size > 1
        set_zonal_boundary_rays!(state)
    end
    if y_size > 1
        set_meridional_boundary_rays!(state)
    end

    return
end

function set_boundary_rays!(
    state::State,
    wkb_mode::Union{SingleColumn, MultiColumn},
)
    (; x_size, y_size) = state.namelists.domain

    if x_size > 1
        set_zonal_boundary_rays!(state)
    end
    if y_size > 1
        set_meridional_boundary_rays!(state)
    end
    set_vertical_boundary_rays!(state)

    return
end
