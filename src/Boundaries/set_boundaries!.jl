"""
```julia
set_boundaries!(state::State, variables::AbstractBoundaryVariables)
```

Enforce all boundary conditions for non-turbulence fields.

```julia 
set_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    turbulence::TKE,
)
```

Enforce all boundary conditions for turbulence fields.

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

function set_boundaries!(state::State, variables::AbstractBoundaryVariables)
    set_zonal_boundaries!(state, variables)
    set_meridional_boundaries!(state, variables)
    set_vertical_boundaries!(state, variables)

    set_tracer_zonal_boundaries!(state, variables)
    set_tracer_meridional_boundaries!(state, variables)
    set_tracer_vertical_boundaries!(state, variables)

    return
end

function set_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    turbulence::TKE,
)
    set_turbulence_zonal_boundaries!(state, variables)
    set_turbulence_meridional_boundaries!(state, variables)
    set_turbulence_vertical_boundaries!(state, variables)

    return
end
