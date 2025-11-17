"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
```

Enforce meridional boundary conditions for tracers by dispatching to a tracer-configuration-specific method.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::AbstractTracer,
)
```

Enforce meridional boundary conditions for tracers.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::AbstractTracer,
)
```

Enforce meridional boundary conditions for reconstructions of tracers.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_model::AbstractWKBMode,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    tracer_setup::AbstractTracer,
)
```

Enforce meridional boundary conditions for tracer-gravity-wave-integral fields.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracer_setup::AbstractTracer,
)
```

Enforce meridional boundary conditions for tracer-gravity-wave-tendency fields.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function set_tracer_meridional_boundaries! end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; tracer_setup) = state.namelists.tracer
    set_tracer_meridional_boundaries!(state, variables, tracer_setup)
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::NoTracer,
)
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::AbstractTracer,
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

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::NoTracer,
)
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::AbstractTracer,
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

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_model::AbstractWKBMode,
    tracer_setup::NoTracer,
)
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
    tracer_setup::AbstractTracer,
)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:uchi, :vchi, :wchi)
        set_meridional_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracer_setup::NoTracer,
)
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
    tracer_setup::AbstractTracer,
)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:dchidt,)
        set_meridional_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
