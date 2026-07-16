"""
```julia
compute_gw_tracer_tendencies!(state::State, i::Integer, j::Integer, k::Integer)
```

Compute the leading-order tracer forcing by dispatching to the appropriate method.

```julia
compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::Val{:NoTracer},
)
```

Return for configurations without tracer transport.

```julia
compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::Val{:TracerOn},
)
```

Compute and return the leading-order tracer forcing at ``\\left(i, j, k\\right)``.

Calculates the tendency that is to be added to the tracer equations, given by

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} \\chi_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{\\bar{\\rho}}\\left[\\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle\\right)_{i + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle\\right)_{i - 1}}{2 \\Delta \\hat{x}} + G^{13} \\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle\\right)_{k + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\hat{z}}\\right.\\\\
    & \\qquad \\qquad + \\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle\\right)_{j + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle\\right)_{j - 1}}{2 \\Delta \\hat{y}} + G^{23} \\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle\\right)_{k + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle\\right)_{k - 1}}{2 \\Delta \\hat{z}}\\\\
    & \\qquad \\qquad + \\left.\\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle\\right)_{k + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle\\right)_{k - 1}}{2 J \\Delta \\hat{z}}\\right] ,
\\end{align*}
```

where ``\\chi_\\mathrm{b}`` is the resolved tracer and ``\\rho_\\mathrm{b}`` is the resolved density (including the reference part ``\\bar{\\rho}``). For a documentation of the fluxes, see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_tracer_integrals!`](@ref).

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `tracer_setup`: General tracer-transport configuration.
"""
function compute_gw_turbulent_tendencies! end

function compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulent_tendencies!(
        state,
        i,
        j,
        k,
        Val(turbulence_scheme),
    )
    return
end

function compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

@ivy function compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
    (; gw_shear) = state.turbulence.turbulencewkbtendencies

    gw_shear[i, j, k] = turbulence_diffusion_coefficient(state, i, j, k, KM()) * 
        gw_shear[i, j, k]

    return
end
