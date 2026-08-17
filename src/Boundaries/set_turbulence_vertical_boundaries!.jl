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

Enforce vertical boundary conditions for turbulence WKB variables by dispatching to the appropriate method.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBIntegrals,
    turbulence_scheme::Val{:TKEScheme},
)
```

Return for turbulence WKB integrals as these are only accessed within the grid cells. (See also [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_turbulence_tendencies!`](@ref))

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Val{:TKEScheme},
)
```

Enforce vertical boundary conditions for turbulence WKB tendencies required for smoothing.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence parameterization configuration.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_turbulence_tendencies!`](@ref)
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
            getfield(turbulencefluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nz == z_size
        for field in fieldnames(TurbulenceFluxes)
            getfield(turbulencefluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme set_turbulence_vertical_boundaries!(
        state::State,
        variables::AbstractBoundaryWKBVariables,
        Val(turbulence_scheme),
    )

    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    turbulence_scheme::Val{:TKEScheme},
)
    return
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Val{:TKEScheme},
)
    (; namelists, domain) = state
    (; turbulencewkbtendencies) = state.turbulence

    for field in fieldnames(TurbulenceWKBTendencies)
        set_vertical_boundaries_of_field!(
            getfield(turbulencewkbtendencies, field),
            namelists,
            domain,
            +,
        )
    end

    return
end
