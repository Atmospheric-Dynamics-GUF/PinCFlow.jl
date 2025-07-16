"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
```
"""
function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
    return
end

"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
```
"""
function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_meridional_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
        )
    end

    return
end

"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
```
"""
function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
    return
end

"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
```
"""
function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_meridional_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
        )
    end

    return
end
