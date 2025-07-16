"""
```julia
IceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```
"""
struct IceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phin::A
    phiq::A
    phiqv::A
end

"""
```julia
IceFluxes(namelists::Namelists, domain::Domain)
```
"""
function IceFluxes(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice

    return IceFluxes(domain, icesetup)
end

"""
```julia
IceFluxes(domain::Domain, icesetup::NoIce)
```
"""
function IceFluxes(domain::Domain, icesetup::NoIce)
    phin = zeros(0, 0, 0, 0)
    phiq = zeros(0, 0, 0, 0)
    phiqv = zeros(0, 0, 0, 0)

    return IceFluxes(phin, phiq, phiqv)
end

"""
```julia
IceFluxes(domain::Domain, icesetup::AbstractIce)
```
"""
function IceFluxes(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    phin = zeros(nxx, nyy, nzz, 3)
    phiq = zeros(nxx, nyy, nzz, 3)
    phiqv = zeros(nxx, nyy, nzz, 3)

    return IceFluxes(phin, phiq, phiqv)
end
