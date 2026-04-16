"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
```

Enforce meridional boundary conditions for tracers by dispatching to the appropriate method.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::Val{:NoTracer},
)
```

Return for configurations without tracer transport.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::Val{:TracerOn},
)
```

Enforce meridional boundary conditions for tracer predictands.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::Val{:TracerOn},
)
```

Enforce meridional boundary conditions for tracer reconstructions.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    tracer_setup::Val{:TracerOn},
)
```

Enforce meridional boundary conditions for tracer WKB quantities by dispatching to the appropriate method.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
```

Enforce meridional boundary conditions for tracer WKB integrals.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
```

Enforce meridional boundary conditions for tracer WKB tendencies.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function set_tracer_meridional_boundaries! end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)
    (; tracer_setup) = state.namelists.tracer
    @dispatch_tracer_setup set_tracer_meridional_boundaries!(
        state,
        variables,
        Val(tracer_setup),
    )
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::Val{:NoTracer},
)
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::Val{:TracerOn},
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
    tracer_setup::Val{:TracerOn},
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
    variables::AbstractBoundaryWKBVariables,
    tracer_setup::Val{:TracerOn},
)
    (; wkb_mode) = state.namelists.wkb
    @dispatch_wkb_mode set_tracer_meridional_boundaries!(
        state,
        variables,
        Val(wkb_mode),
    )
    return
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; tracerwkbintegrals) = state.tracer

    for field in fieldnames(TracerWKBIntegrals)
        set_meridional_boundaries_of_field!(
            getfield(tracerwkbintegrals, field),
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
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    (; namelists, domain) = state
    (; tracerwkbtendencies) = state.tracer

    for field in fieldnames(TracerWKBTendencies)
        set_meridional_boundaries_of_field!(
            getfield(tracerwkbtendencies, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
