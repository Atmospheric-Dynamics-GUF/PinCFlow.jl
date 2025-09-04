"""
```julia 
set_tracer_field_zero!(state)
```

Reset the gravity-wave-induced tracer fluxes to zero by dispatching over tracer configurations.

```julia 
set_tracer_field_zero!(state::State, tracersetup::NoTracer)
```

Return for configurations without tracer transport.

```julia 
set_tracer_field_zero!(state::State, tracersetup::AbstractTracer)
```

Set the gravity-wave-induced tracer fluxes to zero.

# Arguments:

  - `state`: Model state.

  - `tracersetup`: General tracer-transport configuration.

"""
function set_tracer_field_zero! end

function set_tracer_field_zero!(state::State)
    (; tracersetup) = state.namelists.tracer

    set_tracer_field_zero!(state, tracersetup)

    return
end

function set_tracer_field_zero!(state::State, tracersetup::NoTracer)
    return
end

function set_tracer_field_zero!(state::State, tracersetup::AbstractTracer)
    (; chiq0) = state.tracer.tracerforcings

    for field in fieldnames(TracerWKBImpact)
        getfield(chiq0, field) .= 0.0
    end
end
