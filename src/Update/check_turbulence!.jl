"""
```julia
check_turbulence!(state::State)::Nothing
```

Enforce a minimum mass-specific turbulence energy value by dispatching to the turbulence scheme specific method.

```julia
check_turbulence!(state::State, turbulence_scheme::Val{:TKEScheme})::Nothing
```

Enforce a minimum mass-specific turbulent kinetic energy set by `state.turbulence.turbulenceconstants.tkemin`.

```julia
check_turbulence!(state::State, turbulence_scheme::Val{:TTEScheme})::Nothing
```

Enforce a minimum mass-specific turbulent kinetic energy and turbulent available potential energy set by `state.turbulence.turbulenceconstants.tkemin`.

# Arguments

  - `state`: Model state.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
function check_turbulence! end

function check_turbulence!(state::State)::Nothing
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme check_turbulence!(state, Val(turbulence_scheme))
    nothing
end

@ivy function check_turbulence!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
)::Nothing
    (; k0, k1, j0, j1, i0, i1) = state.domain
    (; tke) = state.turbulence.turbulencepredictands
    (; tkemin) = state.turbulence.turbulenceconstants
    (; rhobar) = state.atmosphere

    for k in k0:k1, j in j0:j1, i in i0:i1
        tke[i, j, k] = max(tke[i, j, k], tkemin * rhobar[i, j, k])
    end

    nothing
end

@ivy function check_turbulence!(
    state::State,
    turbulence_scheme::Val{:TTEScheme},
)::Nothing
    (; k0, k1, j0, j1, i0, i1) = state.domain
    (; tke, tape) = state.turbulence.turbulencepredictands
    (; tkemin) = state.turbulence.turbulenceconstants
    (; rhobar) = state.atmosphere

    for k in k0:k1, j in j0:j1, i in i0:i1
        tke[i, j, k] = max(tke[i, j, k], tkemin * rhobar[i, j, k])
        tape[i, j, k] = max(tape[i, j, k], tkemin * rhobar[i, j, k])
    end

    nothing
end
