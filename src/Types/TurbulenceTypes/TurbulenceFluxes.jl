"""
```julia
TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```
"""
struct TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phitke::A
    phitte::A
end

"""
```julia
TurbulenceFluxes(namelists::Namelists, domain::Domain)
```
"""
function TurbulenceFluxes(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceFluxes(domain, turbulencesetup)
end

"""
```julia
TurbulenceFluxes(domain::Domain, turbulencesetup::NoTurbulence)
```
"""
function TurbulenceFluxes(domain::Domain, turbulencesetup::NoTurbulence)
    phitke = zeros(0, 0, 0, 0)
    phitte = zeros(0, 0, 0, 0)

    return TurbulenceFluxes(phitke, phitte)
end

"""
```julia
TurbulenceFluxes(domain::Domain, turbulencesetup::AbstractTurbulence)
```
"""
function TurbulenceFluxes(domain::Domain, turbulencesetup::AbstractTurbulence)
    (; nxx, nyy, nzz) = domain

    phitke = zeros(nxx, nyy, nzz, 3)
    phitte = zeros(nxx, nyy, nzz, 3)

    return TurbulenceFluxes(phitke, phitte)
end
