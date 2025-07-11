struct TracerTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dchi::A
end

function TracerTendencies(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer
    return TracerTendencies(domain, tracersetup)
end

function TracerTendencies(domain::Domain, tracersetup::NoTracer)
    dchi = zeros(0, 0, 0)

    return TracerTendencies(dchi)
end

function TracerTendencies(domain::Domain, tracersetup::AbstractTracer)
    (; nxx, nyy, nzz) = domain

    dchi = zeros(nxx, nyy, nzz)

    return TracerTendencies(dchi)
end
