"""
```julia
set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    model::AbstractModel,
)
```

Return in non-compressible modes.

```julia
set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
```

Enforce vertical boundary conditions for mass-weighted potential temperature in compressible mode (line reflection).

```julia
set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
```

Enforce vertical boundary conditions for vertical mass-weighted potential-temperature flux (no flux through boundaries).

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
  - `model`: Dynamic equations.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_compressible_vertical_boundaries! end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    model::AbstractModel,
)
    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; p) = state.variables.predictands

    set_vertical_boundaries_of_field!(p, namelists, domain, zboundaries, +)

    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; phip) = state.variables.fluxes

    if ko == 0
        phip[:, :, k0 - 1, 3] .= 0.0
    end

    if ko + nzz == sizezz
        phip[:, :, k1, 3] .= 0.0
    end

    return
end
