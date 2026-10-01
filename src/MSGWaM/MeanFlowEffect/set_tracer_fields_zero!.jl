"""
```julia
set_tracer_fields_zero!(state)
```

Reset the gravity-wave-induced tracer fluxes and tracer tendencies to zero by dispatching over tracer configurations.

```julia
set_tracer_fields_zero!(state::State, tracer_setup::Val{:NoTracer})
```

Return for configurations without tracer transport.

```julia
set_tracer_fields_zero!(state::State, tracer_setup::Val{:TracerOn})
```

Set the gravity-wave-induced tracer fluxes and tracer tendencies to zero.

# Arguments

  - `state`: Model state.

  - `tracer_setup`: General tracer-transport configuration.
"""
function set_tracer_fields_zero! end

function set_tracer_fields_zero!(state::State)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup set_tracer_fields_zero!(state, Val(tracer_setup))

    return
end

function set_tracer_fields_zero!(state::State, tracer_setup::Val{:NoTracer})
    return
end

function set_tracer_fields_zero!(state::State, tracer_setup::Val{:TracerOn})
    (; tracerwkbtendencies, tracerwkbintegrals) = state.tracer

    for field_name in fieldnames(TracerWKBTendencies)
        field = getfield(tracerwkbtendencies, field_name)
        @share field = 0.0
    end
    for field_name in fieldnames(TracerWKBIntegrals)
        field = getfield(tracerwkbintegrals, field_name)
        @share field = 0.0
    end

    return
end
