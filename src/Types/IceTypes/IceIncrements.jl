"""
```julia
IceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of ice variables.

```julia
IceIncrements(namelists::Namelists, domain::Domain)::IceIncrements
```

Construct an `IceIncrements` instance with dimensions depending on the general ice-physics configuration, by dispatching to the appropriate method.

```julia
IceIncrements(domain::Domain, ice_setup::NoIce)::IceIncrements
```

Construct an `IceIncrements` instance with zero-size arrays for configurations without ice physics.

```julia
IceIncrements(domain::Domain, ice_setup::AbstractIce)::IceIncrements
```

Construct an `IceIncrements` instance with zero-initialized arrays.

# Fields

  - `dn::A`: Runge-Kutta update of the ice-crystal number concentration.

  - `dq::A`: Runge-Kutta update of the ice mixing ratio.

  - `dqv::A`: Runge-Kutta update of the water-vapor mixing ratio.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `ice_setup`: General ice-physics configuration.
"""
struct IceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
    dn::A
    dq::A
    dqv::A
end

function IceIncrements(namelists::Namelists, domain::Domain)::IceIncrements
    (; ice_setup) = namelists.ice
    return IceIncrements(domain, ice_setup)
end

function IceIncrements(domain::Domain, ice_setup::NoIce)::IceIncrements
    dn = zeros(0, 0, 0)
    dq = zeros(0, 0, 0)
    dqv = zeros(0, 0, 0)

    return IceIncrements(dn, dq, dqv)
end

function IceIncrements(domain::Domain, ice_setup::AbstractIce)::IceIncrements
    (; nxx, nyy, nzz) = domain

    dn = zeros(nxx, nyy, nzz)
    dq = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)

    return IceIncrements(dn, dq, dqv)
end
