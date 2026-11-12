"""
```julia
set_turbulence_fields_zero!(state)::Nothing
```

Reset the gravity-wave shear and the turbulence tendencies to zero.

# Arguments

  - `state`: Model state.
"""
function set_turbulence_fields_zero! end

function set_turbulence_fields_zero!(state::State)::Nothing
    (; turbulencewkbtendencies, turbulencewkbintegrals) = state.turbulence

    for field in fieldnames(TurbulenceWKBIntegrals)
        getfield(turbulencewkbintegrals, field) .= 0.0
    end

    for field in fieldnames(TurbulenceWKBTendencies)
        getfield(turbulencewkbtendencies, field) .= 0.0
    end

    nothing
end
