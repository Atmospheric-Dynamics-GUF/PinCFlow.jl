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

Enforce meridional boundary conditions for the turbulence impact of the gravity-wave shear.

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
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme set_turbulence_meridional_boundaries!(
        state::State,
        variables::AbstractBoundaryWKBVariables,
        Val(turbulence_scheme),
    )

    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    turbulence_scheme::Val{:TKEScheme},
)
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Val{:TKEScheme},
)
    (; wkb_mode) = state.namelists.wkb

    @dispatch_wkb_mode set_turbulence_meridional_boundaries!(
        state,
        variables,
        Val(wkb_mode),
    )
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Val{:NoWKB},
)
    return
end

function set_turbulence_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; turbulencewkbtendencies) = state.turbulence

    for field in fieldnames(TurbulenceWKBTendencies)
        set_meridional_boundaries_of_field!(
            getfield(turbulencewkbtendencies, field),
            namelists,
            domain,
        )
    end

    return
end
