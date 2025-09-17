"""
```julia
backup_predictands(
    state::State
)::Tuple{<:Predictands, <:TracerPredictands}
```

Return a tuple with a copy of the predictands and tracer predictands.

# Arguments 

  - `state`: Model state 

"""

function backup_predictands end 

function backup_predictands(
    state::State
)::Tuple{<:Predictands, <:TracerPredictands}
    p0 = deepcopy(state.variables.predictands)
    chi0 = deepcopy(state.tracer.tracerpredictands) 

    return (p0, chi0)
end