"""
```julia
set_tracer_zonal_boundaries!(state::State, variables::AbstractBoundaryVariables)
```

Enforce zonal boundary conditions for tracers by dispatching to the appropriate method.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::TracerOn,
)
```

Enforce zonal boundary conditions for tracer predictands.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::TracerOn,
)
```

Enforce zonal boundary conditions for tracer reconstructions.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractWKBBoundaryVariables,
    tracer_setup::TracerOn,
)
```

Enforce zonal boundary conditions for tracer WKB quantities by dispatching to the appropriate method.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
```

Enforce zonal boundary conditions for tracer WKB integrals.

```julia
set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
```

Enforce zonal boundary conditions for tracer WKB tendencies.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_tracer_zonal_boundaries! end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; tracer_setup) = state.namelists.tracer
    set_tracer_zonal_boundaries!(state, variables, tracer_setup)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::TracerOn,
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
    tracer_setup::TracerOn,
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
    variables::AbstractWKBBoundaryVariables,
    tracer_setup::TracerOn,
)
    (; wkb_mode) = state.namelists.wkb
    set_tracer_zonal_boundaries!(state, variables, wkb_mode)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
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
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
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
