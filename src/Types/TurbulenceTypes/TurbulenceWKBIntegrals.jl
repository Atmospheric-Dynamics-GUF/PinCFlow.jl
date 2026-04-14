"""
```julia
TurbulenceWKBIntegrals
```

Integrals of ray-volume properties for tracer fluxes.
"""
struct TurbulenceWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    shear::A
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBIntegrals
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme return TurbulenceWKBIntegrals(namelists, domain, Val(turbulence_scheme))
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceWKBIntegrals
    return TurbulenceWKBIntegrals(zeros(0, 0, 0))
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceWKBIntegrals
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TurbulenceWKBIntegrals(domain, Val(wkb_mode))
end

function TurbulenceWKBIntegrals(domain::Domain, wkb_mode::Val{:NoWKB})::TurbulenceWKBIntegrals
    return TurbulenceWKBIntegrals(zeros(0, 0, 0))
end

function TurbulenceWKBIntegrals(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBIntegrals
    (; nxx, nyy, nzz) = domain

    return TurbulenceWKBIntegrals(zeros(nxx, nyy, nzz))
end
