"""
```julia
set_tracer_field_zero!(state)
```

Reset the gravity-wave-induced tracer fluxes to zero by dispatching over tracer configurations.

```julia
set_tracer_field_zero!(state::State, tracer_setup::NoTracer)
```

Return for configurations without tracer transport.

```julia
set_tracer_field_zero!(state::State, tracer_setup::TracerOn)
```

Set the gravity-wave-induced tracer fluxes to zero.

# Arguments:

  - `state`: Model state.

  - `tracer_setup`: General tracer-transport configuration.
"""
function set_tracer_field_zero! end

function set_tracer_field_zero!(state::State)
    (; tracer_setup) = state.namelists.tracer

    set_tracer_field_zero!(state, tracer_setup)

    return
end

function set_tracer_field_zero!(state::State, tracer_setup::NoTracer)
    return
end

function set_tracer_field_zero!(state::State, tracer_setup::TracerOn)
    (; chiq0) = state.tracer.tracerforcings

    for field in fieldnames(TracerWKBImpact)
        getfield(chiq0, field) .= 0.0
    end

    return
end
