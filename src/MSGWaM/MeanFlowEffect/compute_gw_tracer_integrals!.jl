"""
```julia
compute_gw_tracer_integrals!(state::State, parameters::NamedTuple)
```

Compute the leading-order gravity-wave-tracer fluxes by dispatching to the appropriate method.

```julia
compute_gw_tracer_integrals!(
    state::State,
    parameters::NamedTuple,
    tracer_setup::Val{:NoTracer},
)
```

Return for configurations without tracer transport.

```julia
compute_gw_tracer_integrals!(
    state::State,
    parameters::NamedTuple,
    tracer_setup::Val{:TracerOn},
)
```

Compute the leading-order gravity-wave tracer fluxes at ``(i, j, k)``.

The zonal, meridional, and vertical fluxes are given by

```math
\\begin{align*}
    \\bar{\\rho}\\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle & = \\frac{\\bar{\\rho}}{2} \\sum_{r, \\lambda,\\mu,\\nu} \\left[F \\Re \\left(u_\\mathrm{w}\\chi^*_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\bar{\\rho}\\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle & = \\frac{\\bar{\\rho}}{2} \\sum_{r,  \\lambda, \\mu, \\nu} \\left[F \\Re \\left(v_\\mathrm{w} \\chi^*_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\bar{\\rho}\\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle & = \\frac{\\bar{\\rho}}{2} \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(w_\\mathrm{w} \\chi^*_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},
\\end{align*}
```

with flux contributions given by 

```math 
\\begin{align*}
\\frac{\\bar{\\rho}}{2} \\Re \\left(u_{\\mathrm{w}, r}\\chi^*_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r - m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r\\right], \\\\
\\frac{\\bar{\\rho}}{2} \\Re \\left(v_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_{r}} \\frac{m_r}{\\left|\\boldsymbol{k}_{r}\\right|^2} \\mathcal{A}_r \\left[m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r - k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r\\right], \\\\
\\frac{\\bar{\\rho}}{2} \\Re \\left(w_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r - l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r\\right].
\\end{align*}
```

# Arguments

  - `state`: Model state.

  - `parameters`: Named tuple containing parameters used for the computations.

  - `tracer_setup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)
"""
function compute_gw_tracer_integrals! end

function compute_gw_tracer_integrals!(state::State, parameters::NamedTuple)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup compute_gw_tracer_integrals!(
        state,
        parameters,
        Val(tracer_setup),
    )
    return
end

function compute_gw_tracer_integrals!(
    state::State,
    parameters::NamedTuple,
    tracer_setup::Val{:NoTracer},
)
    return
end

@ivy function compute_gw_tracer_integrals!(
    state::State,
    parameters::NamedTuple,
    tracer_setup::Val{:TracerOn},
)
    (; uchi0, vchi0, wchi0) = state.tracer.tracerwkbintegrals
    (; leading_order_impact, next_order_impact) = state.namelists.tracer
    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; thetabar, rhobar) = state.atmosphere
    (; uhat, vhat, what, bhat, pihat, chihat) = state.tracer.tracerwkbamplitudes
    (; kr, lr, mr, fc, omir, dens, factor, xr, yr, zr, iray, jray, kray, n2r) =
        parameters
    (; kappa) = state.constants

    wadr = dens * factor

    if (fc != 0.0 && leading_order_impact) || next_order_impact
        dchidx =
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DX())
        dchidy =
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DY())
        dchidz =
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DZ())
    end

    if fc != 0.0 && leading_order_impact
        coeff = fc / omir * mr * wadr / (kr^2 + lr^2 + mr^2)

        uchi0[iray, jray, kray] += coeff * (lr * dchidz - mr * dchidy)
        vchi0[iray, jray, kray] += coeff * (mr * dchidx - kr * dchidz)
        wchi0[iray, jray, kray] += coeff * (kr * dchidy - lr * dchidx)
    end

    if next_order_impact
        bamp = sqrt(
            dens / rhobar[iray, jray, kray] * 2 * n2r^2 * (kr^2 + lr^2) /
            (omir * (kr^2 + lr^2 + mr^2)),
        )
        uamp =
            1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
            (kr * omir + 1im * lr * fc) *
            bamp
        vamp =
            1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
            (lr * omir - 1im * kr * fc) *
            bamp
        wamp = 1im * omir / n2r * bamp
        bhat[iray, jray, kray] += bamp * factor
        uhat[iray, jray, kray] += uamp * factor
        vhat[iray, jray, kray] += vamp * factor
        what[iray, jray, kray] += wamp * factor
        pihat[iray, jray, kray] +=
            1im / mr * (omir^2 - n2r) / n2r / thetabar[iray, jray, kray] *
            kappa *
            bamp *
            factor
        chihat[iray, jray, kray] +=
            -1im / omir *
            (uamp * dchidx + vamp * dchidy + wamp * dchidz) *
            factor
    end

    return
end
