"""
```julia
backup_predictands(state::State)::Tuple{<:Predictands, <:TracerPredictands}
```

Return a tuple with a copy of the predictands and tracer predictands.

# Arguments

  - `state`: Model state
"""
function backup_predictands end

function backup_predictands(
    state::State,
)::Tuple{<:Predictands, <:TracerPredictands, <:IcePredictands}
    p0 = deepcopy(state.variables.predictands)
    chi0 = deepcopy(state.tracer.tracerpredictands)

    ice_data = state.ice.icepredictands
    active_struct = get_IceActiveVars(ice_data)
    active_fields = ice_active_vars_tuple(active_struct)

    all_fields = fieldnames(typeof(ice_data))
    
    # Iterate over all active fields
    copied_args = map(all_fields) do f
        val = getfield(ice_data, f)

        # only copy active ice fields
        if f in active_fields
            return copy(val)
        else
            return val
        end
    end

    # feed the 6 items in 'copied_args' 
    ice0 = typeof(ice_data)(copied_args...)

    return (p0, chi0, ice0)
end
