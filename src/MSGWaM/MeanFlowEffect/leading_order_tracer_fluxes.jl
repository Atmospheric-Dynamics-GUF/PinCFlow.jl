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
\\overline{\\rho} u_{w,\\alpha}\\chi^*_{w,\\alpha} = \\frac{f}{\\widehat{\\omega}_{\\alpha}}\\frac{m_{\\alpha}}{|\\mathbf{k}_{\\alpha}|^2}\\mathcal{A}_{\\alpha}\\left(l_{\\alpha}\\partial_{z}\\chi_b - m_{\\alpha}\\partial_y\\chi_b \\right).
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
\\overline{\\rho} v_{w,\\alpha}\\chi^*_{w,\\alpha} = \\frac{f}{\\widehat{\\omega}_{\\alpha}}\\frac{m_{\\alpha}}{|\\mathbf{k}_{\\alpha}|^2}\\mathcal{A}_{\\alpha}\\left(m_{\\alpha}\\partial_{x}\\chi_b - k_{\\alpha}\\partial_z\\chi_b \\right).
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
\\overline{\\rho} w_{w,\\alpha}\\chi^*_{w,\\alpha} = \\frac{f}{\\widehat{\\omega}_{\\alpha}}\\frac{m_{\\alpha}}{|\\mathbf{k}_{\\alpha}|^2}\\mathcal{A}_{\\alpha}\\left(k_{\\alpha}\\partial_{y}\\chi_b - l_{\\alpha}\\partial_x\\chi_b \\right).
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
