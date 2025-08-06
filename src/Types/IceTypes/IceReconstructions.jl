"""
```julia
IceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

```julia
IceReconstructions(namelists::Namelists, domain::Domain)
```

```julia
IceReconstructions(domain::Domain, icesetup::NoIce)
```

```julia
IceReconstructions(domain::Domain, icesetup::AbstractIce)
```
"""
struct IceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    ntilde::A
    qtilde::A
    qvtilde::A
end

function IceReconstructions(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice

    return IceReconstructions(domain, icesetup)
end

function IceReconstructions(domain::Domain, icesetup::NoIce)
    ntilde = zeros(0, 0, 0, 0, 0)
    qtilde = zeros(0, 0, 0, 0, 0)
    qvtilde = zeros(0, 0, 0, 0, 0)

    return IceReconstructions(ntilde, qtilde, qvtilde)
end

function IceReconstructions(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    ntilde = zeros(nxx, nyy, nzz, 3, 2)
    qtilde = zeros(nxx, nyy, nzz, 3, 2)
    qvtilde = zeros(nxx, nyy, nzz, 3, 2)

    return IceReconstructions(ntilde, qtilde, qvtilde)
end
