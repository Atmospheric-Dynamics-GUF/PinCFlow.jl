"""
```julia
TracerWKBIntegrals
```

Integrals of ray-volume properties for tracer fluxes.
"""
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

    @dispatch_tracer_setup return TracerWKBIntegrals(namelists, domain, Val(tracer_setup))
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:7]...)
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBIntegrals
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBIntegrals(domain, Val(wkb_mode))
end

function TracerWKBIntegrals(domain::Domain, wkb_mode::Val{:NoWKB})::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:7]...)
end

function TracerWKBIntegrals(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
    (; nxx, nyy, nzz) = domain

    return TracerWKBIntegrals([zeros(nxx, nyy, nzz) for i in 1:7]...)
end
