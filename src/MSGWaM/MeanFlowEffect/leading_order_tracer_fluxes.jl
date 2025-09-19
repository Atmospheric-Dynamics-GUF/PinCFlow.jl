"""
```julia
leading_order_tracer_fluxes(
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
    direction::UChi,
)::AbstractFloat
```

Compute and return the contribution of a ray volume located at `(xlc,ylc,zlc)` to the zonal leading-order gravity-wave tracer fluxes ``\\overline{\\rho}\\langle u'\\chi'\\rangle``.

The flux-contributions are given by

```math
\\overline{\\rho} u_{\\mathrm{w}, r}\\chi^*_{\\mathrm{w}, r} = \\frac{f}{\\widehat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r - m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r\\right].
```

```julia
leading_order_tracer_fluxes(
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
    direction::VChi,
)::AbstractFloat
```

Compute and return the contribution of a ray volume located at `(xlc,ylc,zlc)` to the meridional leading-order gravity-wave tracer fluxes ``\\overline{\\rho}\\langle v'\\chi'\\rangle``.

The fluxes are given by

```math
\\overline{\\rho} v_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r} = \\frac{f}{\\widehat{\\omega}_{r}} \\frac{m_r}{\\left|\\boldsymbol{k}_{r}\\right|^2} \\mathcal{A}_r \\left[m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r - k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r\\right].
```

```julia
leading_order_tracer_fluxes(
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
    direction::WChi,
)::AbstractFloat
```

Compute and return the contribution of a ray volume located at `(xlc,ylc,zlc)` to the vertical leading-order gravity-wave tracer fluxes ``\\overline{\\rho}\\langle w'\\chi'\\rangle``.

The fluxes are given by

```math
\\overline{\\rho} w_{\\mathrm{w}, r} \\chi^*_{\\mathrm{w}, r} = \\frac{f}{\\widehat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r - l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r\\right].
```
"""
function leading_order_tracer_fluxes end

function leading_order_tracer_fluxes(
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
    direction::UChi,
)::AbstractFloat
    dchidy = interpolate_mean_flow(xlc, ylc, zlc, state, DChiDY())
    dchidz = interpolate_mean_flow(xlc, ylc, zlc, state, DChiDZ())

    coeff = fc / omir * wnrm * wadr / (wnrk^2.0 + wnrl^2.0 + wnrm^2.0)

    return coeff * (wnrl * dchidz - wnrm * dchidy)
end

function leading_order_tracer_fluxes(
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
    direction::VChi,
)::AbstractFloat
    dchidx = interpolate_mean_flow(xlc, ylc, zlc, state, DChiDX())
    dchidz = interpolate_mean_flow(xlc, ylc, zlc, state, DChiDZ())

    coeff = fc / omir * wnrm * wadr / (wnrk^2.0 + wnrl^2.0 + wnrm^2.0)

    return coeff * (wnrm * dchidx - wnrk * dchidz)
end

function leading_order_tracer_fluxes(
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
    direction::WChi,
)::AbstractFloat
    dchidx = interpolate_mean_flow(xlc, ylc, zlc, state, DChiDX())
    dchidy = interpolate_mean_flow(xlc, ylc, zlc, state, DChiDY())

    coeff = fc / omir * wnrm * wadr / (wnrk^2.0 + wnrl^2.0 + wnrm^2.0)

    return coeff * (wnrk * dchidy - wnrl * dchidx)
end
