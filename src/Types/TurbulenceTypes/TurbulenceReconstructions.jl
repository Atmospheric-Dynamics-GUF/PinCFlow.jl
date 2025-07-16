"""
```julia
TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```
"""
struct TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    tketilde::A
    ttetilde::A
end

"""
```julia
TurbulenceReconstructions(namelists::Namelists, domain::Domain)
```
"""
function TurbulenceReconstructions(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceReconstructions(domain, turbulencesetup)
end

"""
```julia
TurbulenceReconstructions(domain::Domain, turbulencesetup::NoTurbulence)
```
"""
function TurbulenceReconstructions(
    domain::Domain,
    turbulencesetup::NoTurbulence,
)
    tketilde = zeros(0, 0, 0, 0, 0)
    ttetilde = zeros(0, 0, 0, 0, 0)

    return TurbulenceReconstructions(tketilde, ttetilde)
end

"""
```julia
TurbulenceReconstructions(domain::Domain, turbulencesetup::AbstractTurbulence)
```
"""
function TurbulenceReconstructions(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)
    (; nxx, nyy, nzz) = domain

    tketilde = zeros(nxx, nyy, nzz, 3, 2)
    ttetilde = zeros(nxx, nyy, nzz, 3, 2)

    return TurbulenceReconstructions(tketilde, ttetilde)
end
