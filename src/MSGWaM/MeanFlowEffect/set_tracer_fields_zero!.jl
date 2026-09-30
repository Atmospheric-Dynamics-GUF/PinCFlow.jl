"""
```julia
set_tracer_fields_zero!(state)::Nothing
```

Reset the gravity-wave-induced tracer fluxes and tracer tendencies to zero by dispatching over tracer configurations.

```julia
set_tracer_fields_zero!(state::State, tracer_setup::Val{:NoTracer})::Nothing
```

Return for configurations without tracer transport.

```julia
set_tracer_fields_zero!(state::State, tracer_setup::Val{:TracerOn})::Nothing
```

Set the gravity-wave-induced tracer fluxes and tracer tendencies to zero.

# Arguments

  - `state`: Model state.

  - `tracer_setup`: General tracer-transport configuration.
"""
function set_tracer_fields_zero! end

function set_tracer_fields_zero!(state::State)::Nothing
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup set_tracer_fields_zero!(state, Val(tracer_setup))

    nothing
end

function set_tracer_fields_zero!(
    state::State,
    tracer_setup::Val{:NoTracer},
)::Nothing
    nothing
end

function set_tracer_fields_zero!(
    state::State,
    tracer_setup::Val{:TracerOn},
)::Nothing
    (; tracerwkbtendencies, tracerwkbintegrals) = state.tracer

    for field in fieldnames(TracerWKBTendencies)
        getfield(tracerwkbtendencies, field) .= 0.0
    end
    for field in fieldnames(TracerWKBIntegrals)
        getfield(tracerwkbintegrals, field) .= 0.0
    end

    nothing
end
