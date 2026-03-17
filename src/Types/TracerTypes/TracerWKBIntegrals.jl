struct TracerWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    uchi0::A
    vchi0::A
    wchi0::A
    uchi1::A
    vchi1::A
    wchi1::A
    qchi::A
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBIntegrals
    (; tracer_setup) = namelists.tracer

    return TracerWKBIntegrals(namelists, domain, tracer_setup)
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::NoTracer,
)::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:7]...)
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::TracerOn,
)::TracerWKBIntegrals
    (; wkb_mode) = namelists.wkb

    return TracerWKBIntegrals(domain, wkb_mode)
end

function TracerWKBIntegrals(domain::Domain, wkb_mode::NoWKB)::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:7]...)
end

function TracerWKBIntegrals(
    domain::Domain,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::TracerWKBIntegrals
    (; nxx, nyy, nzz) = domain

    return TracerWKBIntegrals([zeros(nxx, nyy, nzz) for i in 1:7]...)
end