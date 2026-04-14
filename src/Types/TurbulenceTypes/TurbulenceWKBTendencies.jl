"""
```julia
TurbulenceWKBTendencies
```

Tracer tendencies due to gravity waves and turbulence.
"""
struct TurbulenceWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dtkedt::A
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBTendencies
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme return TurbulenceWKBTendencies(
        namelists,
        domain,
        Val(turbulence_scheme),
    )
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceWKBTendencies
    return TurbulenceWKBTendencies(zeros(0, 0, 0))
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceWKBTendencies
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TurbulenceWKBTendencies(domain, Val(wkb_mode))
end

function TurbulenceWKBTendencies(
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TurbulenceWKBTendencies
    return TurbulenceWKBTendencies(zeros(0, 0, 0))
end

function TurbulenceWKBTendencies(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBTendencies
    (; nxx, nyy, nzz) = domain

    return TurbulenceWKBTendencies(zeros(nxx, nyy, nzz))
end
