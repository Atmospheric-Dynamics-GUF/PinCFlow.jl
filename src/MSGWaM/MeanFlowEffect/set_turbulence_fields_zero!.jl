"""
```julia
set_turbulence_fields_zero!(state)::Nothing
```

Reset the gravity-wave shear and the turbulence tendencies to zero by dispatching to the appropriate method.

```julia
set_turbulence_fields_zero!(
    state::State,
    turbulence_setup::Val{:Noturbulence},
)::Nothing
```

Return for configurations without turbulence parameterization.

```julia
set_turbulence_fields_zero!(
    state::State,
    turbulence_setup::Val{:turbulenceOn},
)::Nothing
```

Set the gravity-wave shear and turbulence tendencies to zero.

# Arguments

  - `state`: Model state.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
function set_turbulence_fields_zero! end

function set_turbulence_fields_zero!(state::State)::Nothing
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme set_turbulence_fields_zero!(
        state,
        Val(turbulence_scheme),
    )

    nothing
end

function set_turbulence_fields_zero!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
)::Nothing
    nothing
end

function set_turbulence_fields_zero!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
)::Nothing
    (; turbulencewkbtendencies, turbulencewkbintegrals) = state.turbulence

    for field_name in fieldnames(TurbulenceWKBIntegrals)
        field = getfield(turbulencewkbintegrals, field_name)
        @fuse field = 0.0
    end

    for field_name in fieldnames(TurbulenceWKBTendencies)
        field = getfield(turbulencewkbtendencies, field_name)
        @fuse field = 0.0
    end

    nothing
end
