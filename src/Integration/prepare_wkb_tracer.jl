"""
```julia 
prepare_wkb_tracer(
    state::State,
    operation::Union{MassSpecific, MassWeighted},
)::Nothing
```

Convert the mass-weighted into the mass-specific tracer mixing ratio and vice versa by dispatching to the configuration-specific method.

```julia 
prepare_wkb_tracer(
    state::State,
    operation::Union{MassSpecific, MassWeighted},
    tracer_setup::Val{:NoTracer},
)::Nothing
```

Return for configurations without tracer transport.

```julia 
prepare_wkb_tracer(
    state::State,
    operation::MassSpecific,
    tracer_setup::Val{:TracerOn},
)::Nothing
```

Convert the mass-weighted into the mass-specific tracer mixing ratio.

```julia 
prepare_wkb_tracer(
    state::State,
    operation::MassWeighted,
    tracer_setup::Val{:TracerOn},
)::Nothing
```

Convert the mass-specific into the mass-weighted tracer mixing ratio.

# Arguments 

  - `state`: Model state. 

  - `operation`: Compute the mass-specific or mass-weighted tracer mixing ratio.

  - `tracer_setup`: General tracer-transport configuration.
"""
function prepare_wkb_tracer(
    state::State,
    operation::Union{MassSpecific, MassWeighted},
)::Nothing
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup prepare_wkb_tracer(
        state,
        operation,
        Val(tracer_setup),
    )

    nothing
end

function prepare_wkb_tracer(
    state::State,
    operation::Union{MassSpecific, MassWeighted},
    tracer_setup::Val{:NoTracer},
)::Nothing
    nothing
end

@ivy function prepare_wkb_tracer(
    state::State,
    operation::MassSpecific,
    tracer_setup::Val{:TracerOn},
)::Nothing
    (; tracerpredictands) = state.tracer
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere

    for field in fieldnames(TracerPredictands)
        getfield(tracerpredictands, field) ./= (rho .+ rhobar)
    end

    nothing
end

@ivy function prepare_wkb_tracer(
    state::State,
    operation::MassWeighted,
    tracer_setup::Val{:TracerOn},
)::Nothing
    (; tracerpredictands) = state.tracer
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere

    for field in fieldnames(TracerPredictands)
        getfield(tracerpredictands, field) .*= (rho .+ rhobar)
    end

    nothing
end
