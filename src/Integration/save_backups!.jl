"""
```julia
save_backups!(state::State, variables::Vararg{Symbol})::Nothing
```

Copy the specified fields in `state.variables.predictands` to their counterparts in `state.variables.backups`.

# Arguments

  - `state`: Model state.

  - `variables`: Names of the fields to create backups of.
"""
function save_backups! end

function save_backups!(state::State, variables::Vararg{Symbol})::Nothing
    (; backups, predictands) = state.variables

    for field in variables
        getfield(backups, Symbol(field, :old)) .= getfield(predictands, field)
    end

    nothing
end
