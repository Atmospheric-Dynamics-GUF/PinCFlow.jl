"""
```julia
save_backups!(state::State, variables::Vararg{Symbol})::Nothing
```

Copy the specified fields in `state.variables.predictands`, `state.tracer.tracerpredictands`, or `state.turbulence.turbulencepredictands` to their counterparts in `state.variables.backups`, `state.tracer.tracerbackups`, or `state.turbulence.turbulencebackups`, respectively.

```julia
save_backups!(state::State)::Nothing
```

Copy all fields in `state.variables.predictands`, `state.tracer.tracerpredictands`, and `state.turbulence.turbulencepredictands` to their counterparts in `state.variables.backups`, `state.tracer.tracerbackups`, and `state.turbulence.turbulencebackups`, respectively.

# Arguments

  - `state`: Model state.

  - `variables`: Names of the fields to create backups of.
"""
function save_backups! end

function save_backups!(state::State, variables::Vararg{Symbol})::Nothing
    (; backups, predictands) = state.variables
    (; tracerbackups, tracerpredictands) = state.tracer
    (; turbulencebackups, turbulencepredictands) = state.turbulence

    for field_name in variables
        if field_name in fieldnames(Predictands)
            field = getfield(predictands, field_name)
            backup_field = getfield(backups, Symbol(field_name, :old))
        elseif field_name in fieldnames(TracerPredictands)
            field = getfield(tracerpredictands, field_name)
            backup_field = getfield(tracerbackups, Symbol(field_name, :old))
        elseif field_name in fieldnames(TurbulencePredictands)
            field = getfield(turbulencepredictands, field_name)
            backup_field = getfield(turbulencebackups, Symbol(field_name, :old))
        else
            error("Unknown field!")
        end
        @fuse backup_field = field
    end

    nothing
end

function save_backups!(state::State)::Nothing
    save_backups!(
        state,
        fieldnames(Predictands)...,
        fieldnames(TracerPredictands)...,
        fieldnames(TurbulencePredictands)...,
    )

    nothing
end
