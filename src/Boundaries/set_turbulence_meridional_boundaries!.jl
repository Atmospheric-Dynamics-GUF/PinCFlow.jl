"""
```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
```

Enforce meridional boundary conditions for turbulent kinetic energy.

```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
```

Enforce meridional boundary conditions for reconstructions of turbulent kinetic energy.

```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
```

Enforce meridional boundary conditions for gravity-wave shear.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence parameterization configuration.

# See also

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function set_turbulence_meridional_boundaries! end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
    (; namelists, domain) = state
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_meridional_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_meridional_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
    (; namelists, domain) = state
    (; auxiliaries) = state.wkb

    for field in fieldnames(WKBAuxiliaries)
        set_meridional_boundaries_of_field!(
            getfield(auxiliaries, field),
            namelists,
            domain,
        )
    end

    return
end
