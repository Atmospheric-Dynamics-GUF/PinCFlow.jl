"""
```julia
save_backups!(state::State, variables::Vararg{Symbol})
```

Copy the specified fields in `state.variables.predictands` to their counterparts in `state.variables.backups`.

# Arguments

  - `state`: Model state.

  - `variables`: Names of the fields to create backups of.
"""
function save_backups! end

function save_backups!(state::State, variables::Vararg{Symbol})
    (; backups, predictands) = state.variables

    for field_name in variables
        field = getfield(predictands, field_name)
        backup_field = getfield(backups, Symbol(field_name, :old))
        @share backup_field = field
    end

    return
end
