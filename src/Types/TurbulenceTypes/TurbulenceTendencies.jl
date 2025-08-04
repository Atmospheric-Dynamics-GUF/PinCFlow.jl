"""
```julia
TurbulenceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```
"""
struct TurbulenceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dtke::A
    dtte::A
end

"""
```julia
TurbulenceTendencies(namelists::Namelists, domain::Domain)
```
"""
function TurbulenceTendencies(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceTendencies(domain, turbulencesetup)
end

"""
```julia
TurbulenceTendencies(domain::Domain, turbulencesetup::NoTurbulence)
```
"""
function TurbulenceTendencies(domain::Domain, turbulencesetup::NoTurbulence)
    dtke = zeros(0, 0, 0)
    dtte = zeros(0, 0, 0)

    return TurbulenceTendencies(dtke, dtte)
end

"""
```julia
TurbulenceTendencies(domain::Domain, turbulencesetup::AbstractTurbulence)
```
"""
function TurbulenceTendencies(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)
    (; nxx, nyy, nzz) = domain

    dtke = zeros(nxx, nyy, nzz)
    dtte = zeros(nxx, nyy, nzz)

    return TurbulenceTendencies(dtke, dtte)
end
