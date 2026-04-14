"""
```julia
set_turbulence_field_zero!(state)
```

Reset the gravity-wave-induced tracer fluxes to zero by dispatching over tracer configurations.

```julia
set_turbulence_field_zero!(state::State, tracer_setup::Val{:NoTracer})
```

Return for configurations without tracer transport.

```julia
set_turbulence_field_zero!(state::State, tracer_setup::Val{:TracerOn})
```

Set the gravity-wave-induced tracer fluxes to zero.

# Arguments:

  - `state`: Model state.

  - `turbulence_scheme`: General tracer-transport configuration.
"""
function set_turbulence_field_zero! end

function set_turbulence_field_zero!(state::State)
    (; turbulence_scheme) = state.namelists.tracer

    @dispatch_turbulence_scheme set_turbulence_field_zero!(state, Val(turbulence_scheme))

    return
end

function set_turbulence_field_zero!(state::State, turbulence_scheme::Val{:NoTurbulence})
    return
end

function set_turbulence_field_zero!(state::State, turbulence_scheme::Val{:TKEScheme})
    (; turbulencewkbtendencies, turbulencewkbintegrals) = state.tracer

    for field in fieldnames(TurbulenceWKBTendencies)
        getfield(turbulencewkbtendencies, field) .= 0.0
    end
    for field in fieldnames(TurbulenceWKBIntegrals)
        getfield(turbulencewkbintegrals, field) .= 0.0
    end

    return
end
