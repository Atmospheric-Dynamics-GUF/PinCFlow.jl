"""
```julia
compute_gw_turbulence_integrals!(state::State, parameters::NamedTuple)::Nothing
```

Compute the gravity-wave shear by dispatching to the appropriate method.

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
    parameters::NamedTuple,
)::Nothing
```

Return for configurations without turbulence parameterization.

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    parameters::NamedTuple,
)::Nothing
```

Compute the gravity-wave shear at ``\\left(i, j, k \\right)``, using

```math
    \\mathcal{S}_\\mathrm{w} = \\frac{m_r^4}{\\bar{\\rho}} \\frac{\\left(f^2 + \\hat{\\omega}_r^2\\right)}{\\hat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r.
```

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TTEScheme},
    parameters::NamedTuple,
)::Nothing
```

Compute the gravity-wave shear at ``\\left(i, j, k \\right)`` (see above) and the gravity-wave convective production, using

```math
    \\mathcal{C}_\\mathrm{w} = \\frac{2 m_r^2}{\\bar{\\rho}} \\frac{N_r^2 k_{h, r}^2}{N_r^2 k_{h, r}^2 + f^2m_r^2} \\mathcal{A}_r.
```

# Arguments

  - `state`: Model state.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `parameters`: Wave parameters.
"""
function compute_gw_turbulence_integrals! end

function compute_gw_turbulence_integrals!(
    state::State,
    parameters::NamedTuple,
)::Nothing
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulence_integrals!(
        state,
        Val(turbulence_scheme),
        parameters,
    )

    nothing
end

function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
    parameters::NamedTuple,
)::Nothing
    nothing
end

@ivy function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    parameters::NamedTuple,
)::Nothing
    (; rhobar) = state.atmosphere
    (; gw_shear) = state.turbulence.turbulencewkbintegrals
    (; gw_coupling) = state.namelists.turbulence
    (; fc, wadr, kr, lr, mr, omir, i, j, k) = parameters

    if !gw_coupling
        return nothing
    end

    gw_shear[i, j, k] +=
        mr^4 * wadr / rhobar[i, j, k] * (fc^2 + omir^2) /
        (omir * (kr^2 + lr^2 + mr^2))

    nothing
end

@ivy function compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TTEScheme},
    parameters::NamedTuple,
)::Nothing
    (; rhobar) = state.atmosphere
    (; gw_shear, gw_convective) = state.turbulence.turbulencewkbintegrals
    (; gw_coupling) = state.namelists.turbulence
    (; fc, wadr, n2r, kr, lr, mr, omir, i, j, k) = parameters

    if !gw_coupling
        return nothing
    end

    gw_shear[i, j, k] +=
        mr^4 * wadr / rhobar[i, j, k] * (fc^2 + omir^2) /
        (omir * (kr^2 + lr^2 + mr^2))

    gw_convective[i, j, k] +=
        2 * wadr / omir / rhobar[i, j, k] * mr^2 * n2r * (kr^2 + lr^2) /
        (kr^2 + lr^2 + mr^2)

    nothing
end
