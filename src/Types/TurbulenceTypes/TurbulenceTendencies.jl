"""
```julia
TurbulenceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

```julia
TurbulenceTendencies(namelists::Namelists, domain::Domain)
```

```julia
TurbulenceTendencies(domain::Domain, turbulencesetup::NoTurbulence)
```

```julia
TurbulenceTendencies(domain::Domain, turbulencesetup::AbstractTurbulence)
```
"""
struct TurbulenceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dtke::A
    dtte::A
end

function TurbulenceTendencies(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceTendencies(domain, turbulencesetup)
end

function TurbulenceTendencies(domain::Domain, turbulencesetup::NoTurbulence)
    dtke = zeros(0, 0, 0)
    dtte = zeros(0, 0, 0)

    return TurbulenceTendencies(dtke, dtte)
end

function TurbulenceTendencies(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)
    (; nxx, nyy, nzz) = domain

    dtke = zeros(nxx, nyy, nzz)
    dtte = zeros(nxx, nyy, nzz)

    return TurbulenceTendencies(dtke, dtte)
end
