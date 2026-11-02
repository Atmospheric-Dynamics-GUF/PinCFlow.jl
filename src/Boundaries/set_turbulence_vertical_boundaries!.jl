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

Enforce vertical boundary conditions for turbulence fluxes by dispatching to the vertical boundary condition specific method.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    vertical_boundary_condition::Val{:Periodic},
)::Nothing
```

Return for periodic vertical boundary conditions.

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    vertical_boundary_condition::Val{:SolidWall},
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
    turbulence_scheme::Union{Val{:TKEScheme}, Val{:TTEScheme}},
)::Nothing
```

Return for turbulence WKB integrals as these are not needed in halo/boundary cells (see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_turbulence_tendencies!`](@ref)).

```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Union{Val{:TKEScheme}, Val{:TTEScheme}},
)::Nothing
```

Enforce vertical boundary conditions for turbulence WKB tendencies.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `vertical_boundary_condition`: Vertical boundary conditions.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_turbulence_vertical_boundaries! end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)::Nothing
    (; namelists, domain) = state
    (; tke, tape) = state.turbulence.turbulencepredictands
    (; turbulence_scheme) = state.namelists.turbulence

    set_vertical_boundaries_of_field!(tke, namelists, domain)

    if turbulence_scheme === :TTEScheme
        set_vertical_boundaries_of_field!(tape, namelists, domain)
    end

    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)::Nothing
    (; namelists, domain) = state
    (; tketilde, tapetilde) = state.turbulence.turbulencereconstructions
    (; turbulence_scheme) = state.namelists.turbulence

    set_vertical_boundaries_of_field!(tketilde, namelists, domain)

    if turbulence_scheme === :TTEScheme
        set_vertical_boundaries_of_field!(tapetilde, namelists, domain)
    end

    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
)::Nothing
    (; vertical_boundary_condition) = state.namelists.domain
    @dispatch_vertical_boundary_condition set_turbulence_vertical_boundaries!(
        state,
        variables,
        Val(vertical_boundary_condition),
    )
    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    vertical_boundary_condition::Val{:Periodic},
)::Nothing
    nothing
end

@ivy function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    vertical_boundary_condition::Val{:SolidWall},
)::Nothing
    (; nz, ko, k0, k1) = state.domain
    (; z_size) = state.namelists.domain
    (; phitke, phitape) = state.turbulence.turbulencefluxes
    (; turbulence_scheme) = state.namelists.turbulence

    if ko == 0
        phitke[:, :, k0 - 1, 3] .= 0.0
    end

    if ko + nz == z_size
        phitke[:, :, k1, 3] .= 0.0
    end

    if turbulence_scheme === :TTEScheme
        if ko == 0
            phitape[:, :, k0 - 1, 3] .= 0.0
        end

        if ko + nz == z_size
            phitape[:, :, k1, 3] .= 0.0
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
        state,
        variables,
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
    turbulence_scheme::Union{Val{:TKEScheme}, Val{:TTEScheme}},
)::Nothing
    nothing
end

function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    turbulence_scheme::Union{Val{:TKEScheme}, Val{:TTEScheme}},
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
