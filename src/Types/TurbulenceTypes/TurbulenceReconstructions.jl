"""
```julia
TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

```julia
TurbulenceReconstructions(namelists::Namelists, domain::Domain)
```

```julia
TurbulenceReconstructions(domain::Domain, turbulencesetup::NoTurbulence)
```

```julia
TurbulenceReconstructions(domain::Domain, turbulencesetup::AbstractTurbulence)
```
"""
struct TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    tketilde::A
    ttetilde::A
end

function TurbulenceReconstructions(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceReconstructions(domain, turbulencesetup)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulencesetup::NoTurbulence,
)
    tketilde = zeros(0, 0, 0, 0, 0)
    ttetilde = zeros(0, 0, 0, 0, 0)

    return TurbulenceReconstructions(tketilde, ttetilde)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)
    (; nxx, nyy, nzz) = domain

    tketilde = zeros(nxx, nyy, nzz, 3, 2)
    ttetilde = zeros(nxx, nyy, nzz, 3, 2)

    return TurbulenceReconstructions(tketilde, ttetilde)
end
