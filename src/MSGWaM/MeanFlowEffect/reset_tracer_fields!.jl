"""
```julia
reset_tracer_fields!(state)
```

Reset the gravity-wave-induced tracer fluxes and tracer tendencies to zero by dispatching over tracer configurations.

```julia
reset_tracer_fields!(state::State, tracer_setup::Val{:NoTracer})
```

Return for configurations without tracer transport.

```julia
reset_tracer_fields!(state::State, tracer_setup::Val{:TracerOn})
```

Set the gravity-wave-induced tracer fluxes and tracer tendencies to zero.

# Arguments

  - `state`: Model state.

  - `tracer_setup`: General tracer-transport configuration.
"""
function reset_tracer_fields! end

function reset_tracer_fields!(state::State)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup reset_tracer_fields!(state, Val(tracer_setup))

    return
end

function reset_tracer_fields!(state::State, tracer_setup::Val{:NoTracer})
    return
end

function reset_tracer_fields!(state::State, tracer_setup::Val{:TracerOn})
    (; tracerwkbtendencies, tracerwkbintegrals) =
        state.tracer

    for field in fieldnames(TracerWKBTendencies)
        getfield(tracerwkbtendencies, field) .= 0.0
    end
    for field in fieldnames(TracerWKBIntegrals)
        getfield(tracerwkbintegrals, field) .= 0.0
    end
    
    return
end
