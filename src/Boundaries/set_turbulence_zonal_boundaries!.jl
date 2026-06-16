"""
```julia
set_turbulence_zonal_boundaries!(state::State, variables::BoundaryPredictands)
```

Enforce zonal boundary conditions for turbulent kinetic energy.

```julia
set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
```

Enforce zonal boundary conditions for reconstructions of turbulent kinetic energy.

```julia
set_turbulence_zonal_boundaries!(
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

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_turbulence_zonal_boundaries! end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
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
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:NoTurbulence},
)
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
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; turbulencewkbintegrals) = state.turbulence
    (; wave_impact) = namelists.turbulence

    if wave_impact
        for field in (:gwshear, :gwbuoy)
            set_zonal_boundaries_of_field!(
                getfield(turbulencewkbintegrals, field),
                namelists,
                domain;
                layers = (1, 1, 1),
            )
        end
    end

    return
end

function set_turbulence_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; dtkedt) = state.turbulence.turbulencewkbtendencies
    (; wave_impact) = namelists.turbulence

    if wave_impact
        set_zonal_boundaries_of_field!(
            dtkedt,
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
