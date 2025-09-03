"""
```julia
set_boundaries!(state::State, variables::BoundaryPredictands)
```

Enforce all boundary conditions for predictand fields.

```julia
set_boundaries!(state::State, variables::BoundaryReconstructions)
```

Enforce all boundary conditions for reconstruction fields.

```julia
set_boundaries!(state::State, variables::BoundaryFluxes)
```

Enforce vertical boundary conditions for flux fields (horizontal boundaries are taken care of at the reconstruction stage).

```julia
set_boundaries!(state::State, variables::BoundaryWKBIntegrals)
```

Enforce all boundary conditions for gravity-wave-integral fields.

```julia
set_boundaries!(state::State, variables::BoundaryWKBTendencies)
```

Enforce all boundary conditions for gravity-wave-tendency fields.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_vertical_boundaries!`](@ref)
"""
function set_boundaries! end

function set_boundaries!(state::State, variables::BoundaryPredictands)
    (; zboundaries) = state.namelists.setting

    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

function set_boundaries!(state::State, variables::BoundaryReconstructions)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

function set_boundaries!(state::State, variables::BoundaryFluxes)
    (; zboundaries) = state.namelists.setting
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

function set_boundaries!(state::State, variables::BoundaryWKBIntegrals)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

function set_boundaries!(state::State, variables::BoundaryWKBTendencies)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end
