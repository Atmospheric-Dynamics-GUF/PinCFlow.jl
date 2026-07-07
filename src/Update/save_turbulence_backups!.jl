"""
```julia
save_backups!(state::State, variables::Vararg{Symbol})
```

Copy the specified fields in `state.variables.predictands` to their counterparts in `state.variables.backups`.

# Arguments

  - `state`: Model state.

  - `variables`: Names of the fields to create backups of.
"""
function save_turbulence_backups! end

function save_turbulence_backups!(state::State)
    (; tke) = state.turbulence.turbulencepredictands
    (; tkeold) = state.turbulence.turbulencebackups

    tkeold .= tke
    
    return
end
