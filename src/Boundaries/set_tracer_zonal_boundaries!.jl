"""
```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
```

Enforce zonal boundary conditions for tracers.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
```

Enforce zonal boundary conditions for reconstructions of tracers.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracersetup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_tracer_zonal_boundaries! end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_zonal_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_zonal_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
        )
    end

    return
end
