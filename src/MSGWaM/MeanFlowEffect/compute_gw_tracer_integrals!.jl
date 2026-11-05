"""
```julia
compute_gw_tracer_integrals!(state::State, parameters::NamedTuple)::Nothing
```

Compute the leading-order gravity-wave-tracer fluxes by dispatching to the appropriate method.

```julia
compute_gw_tracer_integrals!(
    state::State,
    tracer_setup::Val{:NoTracer},
    parameters::NamedTuple,
)::Nothing
```

Return for configurations without tracer transport.

```julia
compute_gw_tracer_integrals!(
    state::State,
    tracer_setup::Val{:TracerOn},
    parameters::NamedTuple,
)::Nothing
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
\\frac{\\bar{\\rho}}{2} \\Re \\left(u_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r - m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r\\right], \\\\
\\frac{\\bar{\\rho}}{2} \\Re \\left(v_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_{r}} \\frac{m_r}{\\left|\\boldsymbol{k}_{r}\\right|^2} \\mathcal{A}_r \\left[m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r - k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r\\right], \\\\
\\frac{\\bar{\\rho}}{2} \\Re \\left(w_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r - l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r\\right].
\\end{align*}
```

# Arguments

  - `state`: Model state.

  - `tracer_setup`: General tracer-transport configuration.

  - `parameters`: Wave parameters.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)
"""
function compute_gw_tracer_integrals! end

function compute_gw_tracer_integrals!(
    state::State,
    parameters::NamedTuple,
)::Nothing
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup compute_gw_tracer_integrals!(
        state,
        Val(tracer_setup),
        parameters,
    )
    nothing
end

function compute_gw_tracer_integrals!(
    state::State,
    tracer_setup::Val{:NoTracer},
    parameters::NamedTuple,
)::Nothing
    nothing
end

@ivy function compute_gw_tracer_integrals!(
    state::State,
    tracer_setup::Val{:TracerOn},
    parameters::NamedTuple,
)::Nothing
    (; uchi0, vchi0, wchi0) = state.tracer.tracerwkbintegrals
    (; leading_order_impact) = state.namelists.tracer
    (; i, j, k, fc, kr, lr, mr, wadr, xr, yr, zr, omir) = parameters

    if fc == 0.0 || !leading_order_impact
        return nothing
    end

    coeff = fc / omir * mr * wadr / (kr^2.0 + lr^2.0 + mr^2.0)

    dchidx = interpolate_mean_flow(xr, yr, zr, state, DChiDX())
    dchidy = interpolate_mean_flow(xr, yr, zr, state, DChiDY())
    dchidz = interpolate_mean_flow(xr, yr, zr, state, DChiDZ())

    uchi0[i, j, k] += coeff * (lr * dchidz - mr * dchidy)
    vchi0[i, j, k] += coeff * (mr * dchidx - kr * dchidz)
    wchi0[i, j, k] += coeff * (kr * dchidy - lr * dchidx)

    nothing
end
