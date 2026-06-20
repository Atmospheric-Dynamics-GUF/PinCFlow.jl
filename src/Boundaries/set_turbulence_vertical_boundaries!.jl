"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
```

Enforce vertical boundary conditions for turbulent kinetic energy.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
```

Enforce vertical boundary conditions for reconstructions of turbulent kinetic energy.

```julia
set_turbulence_vertical_boundaries!(state::State, variables::BoundaryFluxes)
```

Set the vertical turbulent kinetic energy fluxes at the vertical boundaries to zero.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
```

Return for WKB-variables.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence parameterization configuration.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_turbulence_vertical_boundaries! end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
    (; namelists, domain) = state
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_vertical_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
            +,
        )
    end

    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

@ivy function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
)
    (; nz, ko, k0, k1) = state.domain
    (; z_size) = state.namelists.domain
    (; turbulencefluxes) = state.turbulence

    if ko == 0
        for field in fieldnames(TurbulenceFluxes)
            @share $getfield(turbulencefluxes, field)[:, :, k0 - 1, 3] = 0.0
        end
    end

    if ko + nz == z_size
        for field in fieldnames(TurbulenceFluxes)
            @share getfield(turbulencefluxes, field)[:, :, k1, 3] = 0.0
        end
    end

    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
    return
end
