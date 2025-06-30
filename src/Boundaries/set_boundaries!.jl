"""
```julia
set_boundaries!(state::State, variables::BoundaryPredictands)
```

Enforce all boundary conditions for predictand fields.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_boundaries!(state::State, variables::BoundaryPredictands)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
```julia
set_boundaries!(state::State, variables::BoundaryReconstructions)
```

Enforce all boundary conditions for reconstruction fields.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_boundaries!(state::State, variables::BoundaryReconstructions)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
```julia
set_boundaries!(state::State, variables::BoundaryFluxes)
```

Enforce vertical boundary conditions for flux fields (horizontal boundaries are taken care of at the reconstruction stage).

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_boundaries!(state::State, variables::BoundaryFluxes)
    (; zboundaries) = state.namelists.setting
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
```julia
set_boundaries!(state::State, variables::BoundaryGWIntegrals)
```

Enforce all boundary conditions for gravity-wave-integral fields.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_boundaries!(state::State, variables::BoundaryGWIntegrals)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end

"""
```julia
set_boundaries!(state::State, variables::BoundaryGWTendencies)
```

Enforce all boundary conditions for gravity-wave-tendency fields.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_boundaries!(state::State, variables::BoundaryGWTendencies)
    (; zboundaries) = state.namelists.setting
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables, zboundaries)
    return
end
