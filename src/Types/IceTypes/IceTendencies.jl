"""
```julia
IceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of ice variables.

```julia
IceTendencies(namelists::Namelists, domain::Domain)
```

Construct an `IceTendencies` instance with dimensions depending on the general ice-physics configuration, by dispatching to the appropriate method.

```julia
IceTendencies(domain::Domain, icesetup::NoIce)
```

Construct an `IceTendencies` instance with zero-size arrays for configurations without ice physics.

```julia
IceTendencies(domain::Domain, icesetup::AbstractIce)
```

Construct an `IceTendencies` instance with zero-initialized arrays.

# Fields

  - `dn::A`: Runge-Kutta update of the ice-crystal number concentration.

  - `dq::A`: Runge-Kutta update of the ice mixing ratio.

  - `dqv::A`: Runge-Kutta update of the water-vapor mixing ratio.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `icesetup`: General ice-physics configuration.
"""
struct IceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dn::A
    dq::A
    dqv::A
end

function IceTendencies(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice
    return IceTendencies(domain, icesetup)
end

function IceTendencies(domain::Domain, icesetup::NoIce)
    dn = zeros(0, 0, 0)
    dq = zeros(0, 0, 0)
    dqv = zeros(0, 0, 0)

    return IceTendencies(dn, dq, dqv)
end

function IceTendencies(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    dn = zeros(nxx, nyy, nzz)
    dq = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)

    return IceTendencies(dn, dq, dqv)
end
