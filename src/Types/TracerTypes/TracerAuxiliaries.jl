struct TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    tracerdiffusion::A
end

function TracerAuxiliaries(
    namelists::Namelists,
    domain::Domain,
)::TracerAuxiliaries
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerAuxiliaries(domain, Val(tracer_setup))
end

function TracerAuxiliaries(
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerAuxiliaries
    return TracerAuxiliaries(
        [zeros(0, 0, 0) for field in fieldnames(TracerAuxiliaries)]...,
    )
end

function TracerAuxiliaries(
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerAuxiliaries
    (; nxx, nyy, nzz) = domain

    return TracerAuxiliaries(
        [zeros(nxx, nyy, nzz) for field in fieldnames(TracerAuxiliaries)]...,
    )
end
