"""
```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
```

Enforce meridional boundary conditions for turbulence energies by dispatching to a turbulence-configuration-specific method.

```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence transport.

```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::TKEScheme,
)
```

Enforce meridional boundary conditions for turbulence energies.

```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::NoTurbulence,
)
```

Return for configurations without turbulence transport.

```julia
set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::TKEScheme,
)
```

Enforce meridional boundary conditions for reconstructions of turbulence energies.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function set_turbulence_meridional_boundaries! end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; turbulence_scheme) = state.namelists.turbulence
    set_turbulence_meridional_boundaries!(state, variables, turbulence_scheme)
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::NoTurbulence,
)
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::TKEScheme,
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
    turbulence_scheme::NoTurbulence,
)
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::TKEScheme,
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
    turbulence_scheme::Union{NoTurbulence, TKEScheme},
)
    return
end
