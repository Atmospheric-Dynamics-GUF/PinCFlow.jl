struct TracerWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dchidt0::A
    dchidt1::A
    dchidtq::A
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBTendencies
    (; tracer_setup) = namelists.tracer

    return TracerWKBTendencies(namelists, domain, tracer_setup)
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::NoTracer,
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::TracerOn,
)::TracerWKBTendencies
    (; wkb_mode) = namelists.wkb

    return TracerWKBTendencies(domain, wkb_mode)
end

function TracerWKBTendencies(
    domain::Domain,
    wkb_mode::NoWKB,
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBTendencies(
    domain::Domain,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::TracerWKBTendencies
    (; nxx, nyy, nzz) = domain

    return TracerWKBTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
