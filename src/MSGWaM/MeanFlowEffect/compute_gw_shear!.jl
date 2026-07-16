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
    \\bar{\\rho}\\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle & = \\frac{\\bar{\\rho}}{2} \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(w_\\mathrm{w} \\chi^*_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu}.
\\end{align*}
```

# Arguments

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
function compute_gw_shear! end

function compute_gw_shear!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_shear!(
        state::State,
        fc::AbstractFloat,
        omir::AbstractFloat,
        wnrk::AbstractFloat,
        wnrl::AbstractFloat,
        wnrm::AbstractFloat,
        wadr::AbstractFloat,
        i::Integer,
        j::Integer,
        k::Integer,
        Val(turbulence_scheme),
    )

    return
end

function compute_gw_shear!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

@ivy function compute_gw_shear!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    wadr::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
    (; rhobar, thetabar) = state.atmosphere
    (; gw_shear) = state.turbulence.turbulencewkbtendencies

    gw_shear[i, j, k] +=
        mr^4  * 
        wadr / 
        rhobar[i, j, k] * 
        (fc^2 + omir^2) / (
            omir * 
            (kr^2 + lr^2 + mr^2)
        )
        
    return
end