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
    ix::Integer,
    jy::Integer,
    kz::Integer,
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
    ix::Integer,
    jy::Integer,
    kz::Integer,
)
    (; uchi, vchi, wchi) = state.tracer.tracerforcings.chiq0

    if fc == 0.0
        return
    end

    uchi[ix, jy, kz] += leading_order_tracer_fluxes(
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
        UCHI(),
    )

    vchi[ix, jy, kz] += leading_order_tracer_fluxes(
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
        VCHI(),
    )

    wchi[ix, jy, kz] += leading_order_tracer_fluxes(
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
        WCHI(),
    )
    return
end