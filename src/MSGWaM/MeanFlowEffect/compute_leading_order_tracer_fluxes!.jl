"""
```julia
compute_leading_order_tracer_fluxes!(
    state::State,
    tracer_setup::NoTracer,
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

Compute the leading-order gravity-wave tracer fluxes by dispatching to the tracer-setup-specific method.

```julia
compute_leading_order_tracer_fluxes!(
    state::State,
    tracer_setup::NoTracer,
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
    tracer_setup::TracerOn,
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
    \\bar{\\rho}\\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle & = \\frac{\\bar{\\rho}}{2} \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(w_\\mathrm{w} \\chi^*_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu}.
\\end{align*}
```

# Arguments:

  - `state`: Model state.

  - `tracer_setup`:  General tracer-transport configuration.

  - `fc`: Coriolis parameter.

  - `omir`: Gravity-wave intrinsic frequency.

  - `wnrk`: Zonal wavenumber.

  - `wnrl`: Meridional wavenumber.

  - `wnrm`: Vertical wavenumber.

  - `wadr`: Phase-space wave-action density.

  - `xlc`: Zonal location of the ray-volume.

  - `ylc`: Meridional location of the ray-volume.

  - `zlc`: Vertical location of the ray-volume.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

# See also:

  - [`PinCFlow.MSGWaM.MeanFlowEffect.leading_order_tracer_fluxes`](@ref)

"""
function compute_leading_order_tracer_fluxes! end

function compute_leading_order_tracer_fluxes!(
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

    compute_leading_order_tracer_fluxes!(
        state,
        tracer_setup,
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
    )
    return
end

function compute_leading_order_tracer_fluxes!(
    state::State,
    tracer_setup::NoTracer,
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
    tracer_setup::TracerOn,
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
    (; uchi0, vchi0, wchi0) = state.tracer.tracerwkbintegrals

    if fc == 0.0
        return
    end

    @ivy uchi0[i, j, k] += leading_order_tracer_fluxes(
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

    @ivy vchi0[i, j, k] += leading_order_tracer_fluxes(
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

    @ivy wchi0[i, j, k] += leading_order_tracer_fluxes(
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
