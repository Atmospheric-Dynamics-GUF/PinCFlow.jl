"""
```julia
set_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)::Nothing
```

Enforce all boundary conditions for the prognostic variables and their reconstructions.

```julia
set_boundaries!(state::State, variables::AbstractBoundaryWKBVariables)::Nothing
```

Enforce all boundary conditions for the gravity wave integral and tendency fields.

```julia
set_boundaries!(state::State, variables::BoundaryFluxes)::Nothing
```

Enforce boundary conditions for flux fields (horizontal boundaries are taken care of at the reconstruction stage).

```julia
set_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
    turbulence::TKE,
)::Nothing
```

Enforce all boundary conditions for turbulence non-flux fields.

```julia
set_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulence::TKE,
)::Nothing
```

Enforce boundary conditions for turbulence flux fields (horizontal boundaries are taken care of at the reconstruction stage).

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
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
)::Nothing
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables)

    set_tracer_zonal_boundaries!(state, variables)
    set_tracer_meridional_boundaries!(state, variables)
    set_tracer_vertical_boundaries!(state, variables)

    nothing
end

function set_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)::Nothing
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables)

    set_tracer_zonal_boundaries!(state, variables)
    set_tracer_meridional_boundaries!(state, variables)
    set_tracer_vertical_boundaries!(state, variables)

    set_turbulence_zonal_boundaries!(state, variables)
    set_turbulence_meridional_boundaries!(state, variables)
    set_turbulence_vertical_boundaries!(state, variables)

    nothing
end

function set_boundaries!(state::State, variables::BoundaryFluxes)::Nothing
    set_vertical_boundaries!(state, variables)

    set_tracer_vertical_boundaries!(state, variables)

    nothing
end

function set_boundaries!(
    state::State,
    variables::Union{BoundaryPredictands, BoundaryReconstructions},
    turbulence::TKE,
)::Nothing
    set_turbulence_zonal_boundaries!(state, variables)
    set_turbulence_meridional_boundaries!(state, variables)
    set_turbulence_vertical_boundaries!(state, variables)

    nothing
end

function set_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulence::TKE,
)::Nothing
    set_turbulence_vertical_boundaries!(state, variables)

    nothing
end
