"""
```julia 
TracerWKBAmplitudes{A <: AbstractArray{<:Complex, 3}}
```

Amplitudes of the unresolved gravity waves.
"""
struct TracerWKBAmplitudes{A <: AbstractArray{<:ComplexF64, 3},}
    uhat::A 
    vhat::A
    what::A
    bhat::A
    pihat::A
    chihat::A
end

function TracerWKBAmplitudes(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBAmplitudes
    (; tracer_setup) = namelists.tracer 

    @dispatch_tracer_setup return TracerWKBAmplitudes(
        namelists, 
        domain,
        Val(tracer_setup),
    )
end

function TracerWKBAmplitudes(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBAmplitudes
    return TracerWKBAmplitudes([zeros(ComplexF64, 0, 0, 0) for i in 1:6]...)
end

function TracerWKBAmplitudes(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBAmplitudes
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBAmplitudes(
        namelists,
        domain,
        Val(wkb_mode),
    )
end

function TracerWKBAmplitudes(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBAmplitudes
    return TracerWKBAmplitudes([zeros(ComplexF64,0, 0, 0) for i in 1:6]...)
end

function TracerWKBAmplitudes(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBAmplitudes
    (; nxx, nyy, nzz) = domain
    (; next_order_impact) = namelists.tracer

    if next_order_impact
        return TracerWKBAmplitudes([zeros(ComplexF64, nxx, nyy, nzz) for i in 1:6]...)
    else 
        return TracerWKBAmplitudes([zeros(ComplexF64,0, 0, 0) for i in 1:6]...)
    end

end
