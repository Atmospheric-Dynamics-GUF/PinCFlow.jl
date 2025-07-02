struct TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    chitilde::A
end

function TracerReconstructions(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer

    return TracerReconstructions(domain, tracersetup)
end

function TracerReconstructions(domain::Domain, tracersetup::NoTracer)
    chitilde = zeros(0, 0, 0, 0, 0)

    return TracerReconstructions(chitilde)
end

function TracerReconstructions(domain::Domain, tracersetup::AbstractTracer)
    (; nxx, nyy, nzz) = domain

    chitilde = zeros(nxx, nyy, nzz, 3, 2)

    return TracerReconstructions(chitilde)
end
