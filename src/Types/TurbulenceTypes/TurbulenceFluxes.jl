"""
```julia
TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

```julia
TurbulenceFluxes(namelists::Namelists, domain::Domain)
```

```julia
TurbulenceFluxes(domain::Domain, turbulencesetup::NoTurbulence)
```

```julia
TurbulenceFluxes(domain::Domain, turbulencesetup::AbstractTurbulence)
```
"""
struct TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phitke::A
    phitte::A
end

function TurbulenceFluxes(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceFluxes(domain, turbulencesetup)
end

function TurbulenceFluxes(domain::Domain, turbulencesetup::NoTurbulence)
    phitke = zeros(0, 0, 0, 0)
    phitte = zeros(0, 0, 0, 0)

    return TurbulenceFluxes(phitke, phitte)
end

function TurbulenceFluxes(domain::Domain, turbulencesetup::AbstractTurbulence)
    (; nxx, nyy, nzz) = domain

    phitke = zeros(nxx, nyy, nzz, 3)
    phitte = zeros(nxx, nyy, nzz, 3)

    return TurbulenceFluxes(phitke, phitte)
end
