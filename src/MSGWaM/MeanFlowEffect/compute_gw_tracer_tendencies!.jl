"""
```julia
compute_gw_tracer_tendencies!(state::State, i::Integer, j::Integer, k::Integer)
```

Compute the tracer forcing by dispatching to the appropriate method.

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

Compute and return the leading- and next-order tracer forcing at ``\\left(i, j, k\\right)^{(0)}``.

Calculates the tendency that is to be added to the tracer equations, given by

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} \\chi_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{w} & = \\left(\\frac{\\partial \\rho_\\mathrm{b} \\chi_\\mathrm{b}}{\\partial t}\\right)^{(0)}_\\mathrm{w} + \\left(\\frac{\\partial \\rho_\\mathrm{b} \\chi_\\mathrm{b}}{\\partial t}\\right)^{(1)}_\\mathrm{w}
\\end{align*}
```

with leading-order impact

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} \\chi_\\mathrm{b}}{\\partial t}\\right)^{(0)}_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{\\bar{\\rho}}\\left[\\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{i + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{i - 1}}{2 \\Delta \\hat{x}} + G^{13} \\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{k + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{k - 1}}{2 \\Delta \\hat{z}}\\right.\\\\
    & \\qquad \\qquad + \\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{j + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{j - 1}}{2 \\Delta \\hat{y}} + G^{23} \\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{k + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{k - 1}}{2 \\Delta \\hat{z}}\\\\
    & \\qquad \\qquad + \\left.\\frac{\\left(\\bar{\\rho} \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{k + 1} - \\left(\\bar{\\rho} \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle^{(0)}\\right)_{k - 1}}{2 J \\Delta \\hat{z}}\\right] 
\\end{align*}
```

and next-order impact

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} \\chi_\\mathrm{b}}{\\partial t}\\right)^{(1)}_\\mathrm{w} & = - \\frac{\\rho_\\mathrm{b}}{2\\bar{\\rho}\\bar{\\theta}}\\left[\\frac{\\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{i + 1} - \\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{i - 1}}{2 \\Delta \\hat{x}} + G^{13} \\frac{\\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{k + 1} - \\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{k - 1}}{2 \\Delta \\hat{z}}\\right.\\\\
    & \\qquad \\qquad + \\frac{\\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{j + 1} - \\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{j - 1}}{2 \\Delta \\hat{y}} + G^{23} \\frac{\\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{k + 1} - \\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{k - 1}}{2 \\Delta \\hat{z}}\\\\
    & \\qquad \\qquad + \\left.\\frac{\\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{k + 1} - \\left(\\bar{\\rho}\\bar{\\theta} \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle^{(1)}\\right)_{k - 1}}{2 J \\Delta \\hat{z}}\\right] .
\\end{align*}
```

Here, ``\\chi_\\mathrm{b}`` is the resolved tracer and ``\\rho_\\mathrm{b}`` is the resolved density (including the reference part ``\\bar{\\rho}``). For a documentation of the fluxes, see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_tracer_integrals!`](@ref).

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `tracer_setup`: General tracer-transport configuration.
"""
function compute_gw_tracer_tendencies! end

function compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup compute_gw_tracer_tendencies!(
        state,
        i,
        j,
        k,
        Val(tracer_setup),
    )
    return
end

function compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::Val{:NoTracer},
)
    return
end

@ivy function compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::Val{:TracerOn},
)
    (; x_size, y_size) = state.namelists.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; uchi0, vchi0, wchi0, uchi1, vchi1, wchi1) =
        state.tracer.tracerwkbintegrals
    (; dchidt0, dchidt1) = state.tracer.tracerwkbtendencies
    (; rho) = state.variables.predictands
    (; rhobar, thetabar) = state.atmosphere
    (; leading_order_impact, next_order_impact) = state.namelists.tracer

    if leading_order_impact
        dchidt0[i, j, k] = 0.0

        if x_size > 1
            dchiu0 =
                (uchi0[i + 1, j, k] - uchi0[i - 1, j, k]) / (2.0 * dx) +
                met[i, j, k, 1, 3] * (uchi0[i, j, k + 1] - uchi0[i, j, k - 1]) /
                (2.0 * dz)
        else
            dchiu0 = 0.0
        end

        if y_size > 1
            dchiv0 =
                (vchi0[i, j + 1, k] - vchi0[i, j - 1, k]) / (2.0 * dy) +
                met[i, j, k, 2, 3] * (vchi0[i, j, k + 1] - vchi0[i, j, k - 1]) /
                (2.0 * dz)
        else
            dchiv0 = 0.0
        end

        dchiw0 =
            (wchi0[i, j, k + 1] - wchi0[i, j, k - 1]) /
            (2.0 * jac[i, j, k] * dz)

        dchidt0[i, j, k] =
            -(rho[i, j, k] + rhobar[i, j, k]) / rhobar[i, j, k] *
            (dchiu0 + dchiv0 + dchiw0)
    end

    if next_order_impact
        dchidt1[i, j, k] = 0.0

        if x_size > 1
            dchiu1 =
                (
                    rhobar[i + 1, j, k] *
                    thetabar[i + 1, j, k] *
                    uchi1[i + 1, j, k] -
                    rhobar[i - 1, j, k] *
                    thetabar[i - 1, j, k] *
                    uchi1[i - 1, j, k]
                ) / (2.0 * dx) +
                met[i, j, k, 1, 3] * (
                    rhobar[i, j, k + 1] *
                    thetabar[i, j, k + 1] *
                    uchi1[i, j, k + 1] -
                    rhobar[i, j, k - 1] *
                    thetabar[i, j, k - 1] *
                    uchi1[i, j, k - 1]
                ) / (2.0 * dz)
        else
            dchiu1 = 0.0
        end

        if y_size > 1
            dchiv1 =
                (
                    rhobar[i, j + 1, k] *
                    thetabar[i, j + 1, k] *
                    vchi1[i, j + 1, k] -
                    rhobar[i, j - 1, k] *
                    thetabar[i, j - 1, k] *
                    vchi1[i, j - 1, k]
                ) / (2.0 * dy) +
                met[i, j, k, 2, 3] * (
                    rhobar[i, j, k + 1] *
                    thetabar[i, j, k + 1] *
                    vchi1[i, j, k + 1] -
                    rhobar[i, j, k - 1] *
                    thetabar[i, j, k - 1] *
                    vchi1[i, j, k - 1]
                ) / (2.0 * dz)
        else
            dchiv1 = 0.0
        end

        dchiw1 =
            (
                rhobar[i, j, k + 1] *
                thetabar[i, j, k + 1] *
                wchi1[i, j, k + 1] -
                rhobar[i, j, k - 1] *
                thetabar[i, j, k - 1] *
                wchi1[i, j, k - 1]
            ) / (2.0 * jac[i, j, k] * dz)

        dchidt1[i, j, k] =
            -(rho[i, j, k] + rhobar[i, j, k]) /
            (2 * rhobar[i, j, k] * thetabar[i, j, k]) *
            (dchiu1 + dchiv1 + dchiw1)
    end

    return
end
