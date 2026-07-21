"""
```julia
set_turbulence_fields_zero!(state)
```

Reset the gravity-wave-induced turbulence fluxes and turbulence tendencies to zero by dispatching over turbulence configurations.

```julia
set_turbulence_fields_zero!(state::State, turbulence_setup::Val{:Noturbulence})
```

Return for configurations without turbulence transport.

```julia
set_turbulence_fields_zero!(state::State, turbulence_setup::Val{:turbulenceOn})
```

Set the gravity-wave-induced turbulence fluxes and turbulence tendencies to zero.

# Arguments

  - `state`: Model state.

  - `turbulence_setup`: General turbulence-transport configuration.
"""
function set_turbulence_fields_zero! end

function set_turbulence_fields_zero!(state::State)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme set_turbulence_fields_zero!(state, Val(turbulence_scheme))

    return
end

function set_turbulence_fields_zero!(
    state::State, 
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

function set_turbulence_fields_zero!(
    state::State, 
    turbulence_scheme::Val{:TKEScheme},
)
    (; turbulencewkbtendencies, turbulencewkbintegrals) = state.turbulence

    for field in fieldnames(TurbulenceWKBIntegrals)
        getfield(turbulencewkbintegrals, field) .= 0.0
    end

    for field in fieldnames(TurbulenceWKBTendencies)
        getfield(turbulencewkbtendencies, field) .= 0.0
    end

    return
end
