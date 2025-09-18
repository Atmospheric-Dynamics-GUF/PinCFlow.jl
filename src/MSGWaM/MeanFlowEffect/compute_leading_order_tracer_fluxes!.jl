"""
```julia
compute_leading_order_tracer_fluxes!(
    state::State,
    tracersetup::NoTracer,
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
compute_leading_order_tracer_fluxes!(
    state::State,
    tracersetup::AbstractTracer,
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
\\overline{\\rho}\\langle u'\\chi'\\rangle & = \\overline{\\rho}\\sum_{\\lambda,\\nu,\\mu,\\alpha}\\left(Fu_{w}\\chi^*_{w}\\right)_{i+\\lambda,j+\\nu,k+\\mu,\\alpha}\\;, \\\\
\\overline{\\rho}\\langle v'\\chi'\\rangle & = \\overline{\\rho}\\sum_{\\lambda,\\nu,\\mu,\\alpha}\\left(Fv_{w}\\chi^*_{w}\\right)_{i+\\lambda,j+\\nu,k+\\mu,\\alpha}\\;, \\\\
\\overline{\\rho}\\langle w'\\chi'\\rangle & = \\overline{\\rho}\\sum_{\\lambda,\\nu,\\mu,\\alpha}\\left(Fw_{w}\\chi^*_{w}\\right)_{i+\\lambda,j+\\nu,k+\\mu,\\alpha}\\;.
\\end{align*}
```

# Arguments:

  - `state`: Model state.

  - `tracersetup`:  General tracer-transport configuration.

  - `fc`: Coriolis parameter.

  - `omir`: Gravity-wave intrinsic frequency.

  - `wnrk`: Zonal wavenumber.

  - `wnrl`: Meridional wavenumber.

  - `wnrm`: Vertical wavenumber.

  - `wadr`: Phase-space wave-action density.

  - `xlc`: Location of the ray-volume in `\\widehat{x}`-direction.

  - `ylc`: Location of the ray-volume in `\\widehat{y}`-direction.

  - `zlc`: Location of the ray-volume in `\\widehat{z}`-direction.

  - `i`: Array index in `\\widehat{x}`-direction.

  - `j`: Array index in `\\widehat{y}`-direction.

  - `k`: Array index in `\\widehat{z}`-direction.

# See also:

  - [`PinCFlow.MSGWaM.MeanFlowEffect.leading_order_tracer_fluxes`](@ref)

"""
function compute_leading_order_tracer_fluxes! end

function compute_leading_order_tracer_fluxes!(
    state::State,
    tracersetup::NoTracer,
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
    return
end

function compute_leading_order_tracer_fluxes!(
    state::State,
    tracersetup::AbstractTracer,
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
    (; uchi, vchi, wchi) = state.tracer.tracerforcings.chiq0

    if fc == 0.0
        return
    end

    @ivy uchi[i, j, k] += leading_order_tracer_fluxes(
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
        UChi(),
    )

    @ivy vchi[i, j, k] += leading_order_tracer_fluxes(
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
        VChi(),
    )

    @ivy wchi[i, j, k] += leading_order_tracer_fluxes(
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
        WChi(),
    )

    return
end
