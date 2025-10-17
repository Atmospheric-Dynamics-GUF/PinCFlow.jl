"""
```julia
set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
```

Enforce vertical boundary conditions for the specified variables in compressible mode by dispatching to a model-specific method.

```julia
set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    model::Union{Boussinesq, PseudoIncompressible},
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

Enforce vertical boundary conditions for the mass-weighted potential temperature in compressible mode (line reflection).

```julia
set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
```

Enforce vertical boundary conditions for the vertical mass-weighted potential-temperature flux (no flux through boundaries).

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
)
    (; model) = state.namelists.atmosphere
    set_compressible_vertical_boundaries!(state, variables, model)
    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    model::Union{Boussinesq, PseudoIncompressible},
)
    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
    (; namelists, domain) = state
    (; p) = state.variables.predictands

    set_vertical_boundaries_of_field!(p, namelists, domain, +)

    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
    (; z_size) = state.namelists.domain
    (; nz, ko, k0, k1) = state.domain
    (; phip) = state.variables.fluxes

    @ivy if ko == 0
        phip[:, :, k0 - 1, 3] .= 0.0
    end

    @ivy if ko + nz == z_size
        phip[:, :, k1, 3] .= 0.0
    end

    return
end
