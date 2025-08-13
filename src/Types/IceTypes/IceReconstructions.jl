"""
```julia
IceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of ice variables.

The first three dimensions represent physical space, the fourth dimension represents the direction in which the reconstruction was performed and the fifth dimension represents the two cell edges of the reconstruction.

```julia
IceReconstructions(namelists::Namelists, domain::Domain)
```

Construct an `IceReconstructions` instance with dimensions depending on the general ice-physics configuration, by dispatching to the appropriate method.

```julia
IceReconstructions(domain::Domain, icesetup::NoIce)
```

Construct an `IceReconstructions` instance with zero-size arrays for configurations without ice physics.

```julia
IceReconstructions(domain::Domain, icesetup::AbstractIce)
```

Construct an `IceReconstructions` instance with zero-initialized arrays.

# Fields

  - `ntilde::A`: Reconstructions of the ice-crystal number concentration.

  - `qtilde::A`: Reconstructions of the ice mixing ratio.

  - `qvtilde::A`: Reconstructions of the water-vapor mixing ratio.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `icesetup`: General ice-physics configuration.
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
