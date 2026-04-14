"""
```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
```

Enforce zonal boundary conditions for turbulence energies by dispatching to a turbulence-configuration-specific method.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::Val{:TKEScheme},
)
```

Enforce zonal boundary conditions for turbulent kinetic energy.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::Val{:TKEScheme},
)
```

Enforce zonal boundary conditions for reconstructions of turbulent kinetic energy.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Union{Val{:NoTurbulence}, Val{:TKEScheme}},
)
```

Return for WKB-variables.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence parameterization configuration.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_turbulence_zonal_boundaries! end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; turbulence_scheme) = state.namelists.turbulence
    @dispatch_turbulence_scheme set_turbulence_zonal_boundaries!(
        state,
        variables,
        Val(turbulence_scheme),
    )
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulence_scheme::Val{:TKEScheme},
)
    (; namelists, domain) = state
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_zonal_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulence_scheme::Val{:TKEScheme},
)
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_zonal_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:TKEScheme},
)
    (; wkb_mode) = state.namelists.wkb
    @dispatch_wkb_mode set_turbulence_zonal_boundaries!(
        state,
        variables,
        Val(wkb_mode),
    )
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Union{Val{:NoTurbulence}, Val{:TKEScheme}},
)
    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; turbulencewkbintegrals) = state.turbulence

    for field in fieldnames(TurbulenceWKBIntegrals)
        set_zonal_boundaries_of_field!(
            getfield(turbulencewkbintegrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; turbulencewkbtendencies) = state.turbulence

    for field in fieldnames(TurbulenceWKBTendencies)
        set_zonal_boundaries_of_field!(
            getfield(turbulencewkbtendencies, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end