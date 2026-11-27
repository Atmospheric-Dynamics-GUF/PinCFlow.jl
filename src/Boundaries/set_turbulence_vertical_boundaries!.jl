"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)::Nothing
```

Enforce vertical boundary conditions for turbulent kinetic energy.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)::Nothing
```

Enforce vertical boundary conditions for reconstructions of turbulent kinetic energy.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
)::Nothing
```

Set the vertical turbulent kinetic energy fluxes at the vertical boundaries to zero.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)::Nothing
```

Enforce vertical boundary conditions for turbulence WKB variables by dispatching to the appropriate method.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:NoTurbulence},
)::Nothing
```

Return for configurations without turbulence parameterization.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBIntegrals,
    turbulence_scheme::Val{:TKEScheme},
)::Nothing
```

Return for turbulence WKB integrals as these are not needed in halo/boundary cells (see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_turbulence_tendencies!`](@ref)).

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Val{:TKEScheme},
)::Nothing
```

Enforce vertical boundary conditions for turbulence WKB tendencies.

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
)::Nothing
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

    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)::Nothing
    (; namelists, domain) = state
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
        )
    end

    nothing
end

@ivy function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
)::Nothing
    (; nz, ko, k0, k1) = state.domain
    (; z_size) = state.namelists.domain
    (; turbulencefluxes) = state.turbulence

    if ko == 0
        for field_name in fieldnames(TurbulenceFluxes)
            field = getfield(turbulencefluxes, field_name)
            @fuse field[:, :, k0 - 1, 3] = 0.0
        end
    end

    if ko + nz == z_size
        for field_name in fieldnames(TurbulenceFluxes)
            field = getfield(turbulencefluxes, field_name)
            @fuse field[:, :, k1, 3] = 0.0
        end
    end

    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)::Nothing
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme set_turbulence_vertical_boundaries!(
        state::State,
        variables::AbstractBoundaryWKBVariables,
        Val(turbulence_scheme),
    )

    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    turbulence_scheme::Val{:NoTurbulence},
)::Nothing
    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    turbulence_scheme::Val{:TKEScheme},
)::Nothing
    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Val{:TKEScheme},
)::Nothing
    (; namelists, domain) = state
    (; gw_coupling) = namelists.turbulence
    (; turbulencewkbtendencies) = state.turbulence

    if !gw_coupling
        return nothing
    end

    for field in fieldnames(TurbulenceWKBTendencies)
        set_vertical_boundaries_of_field!(
            getfield(turbulencewkbtendencies, field),
            namelists,
            domain,
            +,
        )
    end

    nothing
end
