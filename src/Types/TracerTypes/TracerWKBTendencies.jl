"""
```julia
TracerWKBTendencies
```

Tracer tendencies due to gravity waves and turbulence.
"""
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

    @dispatch_tracer_setup return TracerWKBTendencies(
        namelists,
        domain,
        Val(tracer_setup),
    )
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBTendencies
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBTendencies(domain, Val(wkb_mode))
end

function TracerWKBTendencies(
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBTendencies(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBTendencies
    (; nxx, nyy, nzz) = domain

    return TracerWKBTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
