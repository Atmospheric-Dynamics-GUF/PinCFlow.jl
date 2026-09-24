"""
```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
)::Nothing
```

Enforce meridional boundary conditions for tracers by dispatching to the appropriate method.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::Val{:NoTracer},
)::Nothing
```

Return for configurations without tracer transport.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::Val{:TracerOn},
)::Nothing
```

Enforce meridional boundary conditions for tracer predictands.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::Val{:TracerOn},
)::Nothing
```

Enforce meridional boundary conditions for tracer reconstructions.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    tracer_setup::Val{:TracerOn},
)::Nothing
```

Enforce meridional boundary conditions for tracer WKB quantities by dispatching to the appropriate method.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::Nothing
```

Enforce meridional boundary conditions for tracer WKB integrals.

```julia
set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::Nothing
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
)::Nothing
    (; tracer_setup) = state.namelists.tracer
    @dispatch_tracer_setup set_tracer_meridional_boundaries!(
        state,
        variables,
        Val(tracer_setup),
    )
    nothing
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    tracer_setup::Val{:NoTracer},
)::Nothing
    nothing
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracer_setup::Val{:TracerOn},
)::Nothing
    (; namelists, domain) = state
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_meridional_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracer_setup::Val{:TracerOn},
)::Nothing
    (; namelists, domain) = state
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_meridional_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
        )
    end

    nothing
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
    tracer_setup::Val{:TracerOn},
)::Nothing
    (; wkb_mode) = state.namelists.wkb
    @dispatch_wkb_mode set_tracer_meridional_boundaries!(
        state,
        variables,
        Val(wkb_mode),
    )
    nothing
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::Nothing
    (; namelists, domain) = state
    (; tracerwkbintegrals) = state.tracer
    (; leading_order_impact) = namelists.tracer

    if leading_order_impact
        for field in (:uchi0, :vchi0, :wchi0)
            set_meridional_boundaries_of_field!(
                getfield(tracerwkbintegrals, field),
                namelists,
                domain;
                layers = (1, 1, 1),
            )
        end
    end

    nothing
end

function set_tracer_meridional_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::Nothing
    (; namelists, domain) = state
    (; dchidt0) = state.tracer.tracerwkbtendencies
    (; leading_order_impact) = namelists.tracer

    if leading_order_impact
        set_meridional_boundaries_of_field!(
            dchidt0,
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    nothing
end
