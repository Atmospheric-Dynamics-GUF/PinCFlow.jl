"""
```julia
save_backups!(state::State, variables::Vararg{Symbol})
```

Save backup copies of specified predictand variables.

Creates backup copies of the current values of specified predicted variables
by copying them to corresponding backup fields with "old" suffix.

# Arguments

  - `state::State`: Simulation state containing variables and backups
  - `variables::Symbol...`: Variable names to backup (e.g., `:rho`, `:u`, `:v`, `:w`)

# Details

For each specified variable `field`, copies `predictands.field` to `backups.fieldold`.

# Example

```julia
save_backups!(state, :rho, :u, :v, :w)  # Backup density and velocity fields
```
"""
function save_backups!(state::State, variables::Vararg{Symbol})
    (; backups, predictands) = state.variables

    for field in variables
        getfield(backups, Symbol(field, :old)) .= getfield(predictands, field)
    end

    return
end
