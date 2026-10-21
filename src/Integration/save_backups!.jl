"""
```julia
save_backups!(state::State, variables::Vararg{Symbol})
```

Copy the specified fields in `state.variables.predictands`, `state.tracer.tracerpredictands`, or `state.turbulence.turbulencepredictands` to their counterparts in `state.variables.backups`, `state.tracer.tracerbackups`, or `state.turbulence.turbulencebackups`, respectively.

```julia
save_backups!(state::State)
```

Copy all fields in `state.variables.predictands`, `state.tracer.tracerpredictands`, and `state.turbulence.turbulencepredictands` to their counterparts in `state.variables.backups`, `state.tracer.tracerbackups`, and `state.turbulence.turbulencebackups`, respectively.

# Arguments

  - `state`: Model state.

  - `variables`: Names of the fields to create backups of.
"""
function save_backups! end

function save_backups!(state::State, variables::Vararg{Symbol})
    (; backups, predictands) = state.variables
    (; tracerbackups, tracerpredictands) = state.tracer
    (; turbulencebackups, turbulencepredictands) = state.turbulence

    for field in variables
        if field in fieldnames(Predictands)
            getfield(backups, Symbol(field, :old)) .=
                getfield(predictands, field)
        elseif field in fieldnames(TracerPredictands)
            getfield(tracerbackups, Symbol(field, :old)) .=
                getfield(tracerpredictands, field)
        elseif field in fieldnames(TurbulencePredictands)
            getfield(turbulencebackups, Symbol(field, :old)) .=
                getfield(turbulencepredictands, field)
        else
            error("Unknown field!")
        end
    end

    return
end

function save_backups!(state::State)
    save_backups!(
        state,
        fieldnames(Predictands)...,
        fieldnames(TracerPredictands)...,
        fieldnames(TurbulencePredictands)...,
    )

    return
end
