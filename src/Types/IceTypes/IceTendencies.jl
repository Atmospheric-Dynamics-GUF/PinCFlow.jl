"""
```julia
IceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```
"""
struct IceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dn::A
    dq::A
    dqv::A
end

"""
```julia
IceTendencies(namelists::Namelists, domain::Domain)
```
"""
function IceTendencies(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice
    return IceTendencies(domain, icesetup)
end

"""
```julia
IceTendencies(domain::Domain, icesetup::NoIce)
```
"""
function IceTendencies(domain::Domain, icesetup::NoIce)
    dn = zeros(0, 0, 0)
    dq = zeros(0, 0, 0)
    dqv = zeros(0, 0, 0)

    return IceTendencies(dn, dq, dqv)
end

"""
```julia
IceTendencies(domain::Domain, icesetup::AbstractIce)
```
"""
function IceTendencies(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    dn = zeros(nxx, nyy, nzz)
    dq = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)

    return IceTendencies(dn, dq, dqv)
end
