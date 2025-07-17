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
    direction::UCHI,
) 
    dchidy = compute_derivatives(xlc, ylc, zlc, state, DCHIDY())
    dchidz = compute_derivatives(xlc, ylc, zlc, state, DCHIDZ())

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
    direction::VCHI,
)
    dchidx = compute_derivatives(xlc, ylc, zlc, state, DCHIDX())
    dchidz = compute_derivatives(xlc, ylc, zlc, state, DCHIDZ())

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
    direction::WCHI,
)
    dchidx = compute_derivatives(xlc, ylc, zlc, state, DCHIDX())
    dchidy = compute_derivatives(xlc, ylc, zlc, state, DCHIDY())

    coeff = fc / omir * wnrm * wadr / (wnrk^2.0 + wnrl^2.0 + wnrm^2.0)

    return coeff * (wnrk * dchidy - wnrl * dchidx)
end
