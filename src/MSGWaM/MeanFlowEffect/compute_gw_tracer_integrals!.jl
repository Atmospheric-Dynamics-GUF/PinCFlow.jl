"""
```julia
compute_gw_tracer_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Compute the leading-order gravity-wave-tracer fluxes by dispatching to the appropriate method.

```julia
compute_gw_tracer_integrals!(
    state::State,
    tracer_setup::Val{:NoTracer},
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Return for configurations without tracer transport.

```julia
compute_gw_tracer_integrals!(
    state::State,
    tracer_setup::Val{:TracerOn},
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
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

  - `tracer_setup`: General tracer-transport configuration.

  - `fc`: Coriolis parameter.

  - `omir`: Gravity-wave intrinsic frequency.

  - `wnrk`: Zonal wavenumber.

  - `wnrl`: Meridional wavenumber.

  - `wnrm`: Vertical wavenumber.

  - `wadr`: Contributing fraction of the physical-space wave-action density.

  - `xlc`: Zonal location of the ray-volume.

  - `ylc`: Meridional location of the ray-volume.

  - `zlc`: Vertical location of the ray-volume.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)
"""
function compute_gw_tracer_integrals! end

function compute_gw_tracer_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup compute_gw_tracer_integrals!(
        state,
        fc,
        omir,
        wnrk,
        wnrl,
        wnrm,
        wadr,
        xlc,
        ylc,
        zlc,
        i,
        j,
        k,
        Val(tracer_setup),
    )
    return
end

function compute_gw_tracer_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::Val{:NoTracer},
)
    return
end

@ivy function compute_gw_tracer_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::Val{:TracerOn},
)
    (; uchi0, vchi0, wchi0) = state.tracer.tracerwkbintegrals
    (; leading_order_impact) = state.namelists.tracer
    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere

    if fc == 0.0 || !leading_order_impact
        return
    end

    coeff = fc / omir * wnrm * wadr / (wnrk^2.0 + wnrl^2.0 + wnrm^2.0)

    dchidx =
        interpolate_scalar(state, xlc, ylc, zlc, chi ./ (rho .+ rhobar), DX())
    dchidy =
        interpolate_scalar(state, xlc, ylc, zlc, chi ./ (rho .+ rhobar), DY())
    dchidz =
        interpolate_scalar(state, xlc, ylc, zlc, chi ./ (rho .+ rhobar), DZ())

    uchi0[i, j, k] += coeff * (wnrl * dchidz - wnrm * dchidy)
    vchi0[i, j, k] += coeff * (wnrm * dchidx - wnrk * dchidz)
    wchi0[i, j, k] += coeff * (wnrk * dchidy - wnrl * dchidx)

    return
end
