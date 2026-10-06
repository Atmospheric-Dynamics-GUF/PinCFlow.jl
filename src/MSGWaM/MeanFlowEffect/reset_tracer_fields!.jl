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
    (; tracerwkbtendencies, tracerwkbintegrals) = state.tracer

    for field in fieldnames(TracerWKBTendencies)
        getfield(tracerwkbtendencies, field) .= 0.0
    end
    for (field, fieldold) in zip(
        (:uhat, :vhat, :what, :bhat, :chihat),
        (:uhatold, :vhatold, :whatold, :bhatold, :chihatold),
    )
        getfield(tracerwkbintegrals, fieldold) .=
            getfield(tracerwkbintegrals, field)
    end
    for field in (
        :uchi0,
        :vchi0,
        :wchi0,
        :uchi1,
        :vchi1,
        :wchi1,
        :uhat,
        :vhat,
        :what,
        :bhat,
        :pihat,
        :chihat,
    )
        getfield(tracerwkbintegrals, field) .= 0.0
    end

    return
end
