"""
```julia
set_compressible_meridional_boundaries!(state::State, model::AbstractModel)
```

Return in non-compressible modes.

```julia
set_compressible_meridional_boundaries!(state::State, model::Compressible)
```

Enforce meridional boundary conditions for mass-weighted potential temperature in compressible mode.

# Arguments

- `state`: Model state.
- `model`: Dynamic equations.

# See also

- [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function set_compressible_meridional_boundaries! end

function set_compressible_meridional_boundaries!(
    state::State,
    model::AbstractModel,
)
    return
end

function set_compressible_meridional_boundaries!(
    state::State,
    model::Compressible,
)
    (; namelists, domain) = state
    (; p) = state.variables.predictands
    set_meridional_boundaries_of_field!(p, namelists, domain)
    return
end
