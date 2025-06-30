"""
```julia
set_compressible_zonal_boundaries!(state::State, model::AbstractModel)
```

Return in non-compressible modes.

# Arguments

  - `state`: Model state.
  - `model`: Dynamic equations.
"""
function set_compressible_zonal_boundaries!(state::State, model::AbstractModel)
    return
end

"""
```julia
set_compressible_zonal_boundaries!(state, model::Compressible)
```

Enforce zonal boundary conditions for mass-weighted potential temperature in compressible mode.

# Arguments

  - `state`: Model state.
  - `model`: Dynamic equations.
"""
function set_compressible_zonal_boundaries!(state::State, model::Compressible)
    (; namelists, domain) = state
    (; p) = state.variables.predictands
    set_zonal_boundaries_of_field!(p, namelists, domain)
    return
end
