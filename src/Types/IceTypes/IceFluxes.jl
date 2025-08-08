"""
```julia
IceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

Arrays for fluxes of ice variables.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
IceFluxes(namelists::Namelists, domain::Domain)
```

Construct an `IceFluxes` instance with dimensions depending on the general ice-physics configuration, by dispatching to the appropriate method.

```julia
IceFluxes(domain::Domain, icesetup::NoIce)
```

Construct an `IceFluxes` instance with zero-size arrays for configurations without ice physics.

```julia
IceFluxes(domain::Domain, icesetup::AbstractIce)
```

Construct an `IceFluxes` instance with zero-initialized arrays.

# Fields

- `phin::A`: Fluxes of the ice-crystal number concentration.
- `phiq::A`: Fluxes of the ice mixing ratio.
- `phiqv::A`: Fluxes of the water-vapor mixing ratio.

# Arguments

- `namelists`: Namelists with all model parameters.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `icesetup`: General ice-physics configuration.
"""
struct IceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phin::A
    phiq::A
    phiqv::A
end

function IceFluxes(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice

    return IceFluxes(domain, icesetup)
end

function IceFluxes(domain::Domain, icesetup::NoIce)
    phin = zeros(0, 0, 0, 0)
    phiq = zeros(0, 0, 0, 0)
    phiqv = zeros(0, 0, 0, 0)

    return IceFluxes(phin, phiq, phiqv)
end

function IceFluxes(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    phin = zeros(nxx, nyy, nzz, 3)
    phiq = zeros(nxx, nyy, nzz, 3)
    phiqv = zeros(nxx, nyy, nzz, 3)

    return IceFluxes(phin, phiq, phiqv)
end
