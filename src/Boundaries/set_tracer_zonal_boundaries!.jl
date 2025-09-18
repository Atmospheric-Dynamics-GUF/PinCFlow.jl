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

```julia 
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_model::AbstractWKBMode,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia 
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    tracersetup::AbstractTracer,
)
```

Enforce zonal boundary conditions for tracer-gravity-wave-integral fields.

```julia 
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracersetup::NoTracer,
)
```

Return for configurations without tracer transport. 

```julia 
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracersetup::AbstractTracer,
)
```

Enforce zonal boundary conditions for tracer-gravity-wave-tendency fields.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracersetup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MSGWaM.

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

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_model::AbstractWKBMode,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:uchi, :vchi, :wchi)
        set_zonal_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:dchidt,)
        set_zonal_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
