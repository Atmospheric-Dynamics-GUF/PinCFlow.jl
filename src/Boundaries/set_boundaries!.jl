"""
```julia
set_boundaries!(
    state::State,
    variables::Union{
        BoundaryPredictands,
        BoundaryReconstructions,
        AbstractBoundaryWKBVariables,
    },
)
```

Enforce all boundary conditions for non-flux fields.

```julia
set_boundaries!(state::State, variables::BoundaryFluxes)
```

Enforce vertical boundary conditions for flux fields (horizontal boundaries are taken care of at the reconstruction stage).

```julia 
set_boundaries!(
    state::State,
    variables::Union{
        BoundaryPredictands,
        BoundaryReconstructions,
        AbstractBoundaryWKBVariables,
    },
    turbulence::TKE,
)
```

Enforce all boundary conditions for turbulence non-flux fields.

```julia 
set_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulence::TKE,
)
```

Enforce vertical boundary conditions for turbulence flux fields (horizontal boundaries are taken care of at the reconstruction stage).

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `turbulence`: Enforce boundary conditions for turbulence fields.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_vertical_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_tracer_zonal_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_tracer_meridional_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_tracer_vertical_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_turbulence_zonal_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_turbulence_meridional_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_turbulence_vertical_boundaries!`](@ref)
"""
function set_boundaries! end

function set_boundaries!(
    state::State,
    variables::Union{
        BoundaryPredictands,
        BoundaryReconstructions,
        AbstractBoundaryWKBVariables,
    },
)
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables)

    set_tracer_zonal_boundaries!(state, variables)
    set_tracer_meridional_boundaries!(state, variables)
    set_tracer_vertical_boundaries!(state, variables)

    return
end

function set_boundaries!(state::State, variables::BoundaryFluxes)
    set_vertical_boundaries!(state, variables)

    set_tracer_vertical_boundaries!(state, variables)

    return
end

function set_boundaries!(
    state::State,
    variables::Union{
        BoundaryPredictands,
        BoundaryReconstructions,
        AbstractBoundaryWKBVariables,
    },
    turbulence::TKE,
)
    set_turbulence_zonal_boundaries!(state, variables)
    set_turbulence_meridional_boundaries!(state, variables)
    set_turbulence_vertical_boundaries!(state, variables)

    return
end

function set_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulence::TKE,
)
    set_turbulence_vertical_boundaries!(state, variables)

    return
end
