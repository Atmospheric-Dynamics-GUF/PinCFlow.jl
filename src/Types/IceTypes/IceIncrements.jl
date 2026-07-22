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

  - `dN_in`: Runge-Kutta update of the number contration of IN.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `ice_setup`: General ice-physics configuration.
"""
struct IceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
    dn_hom::A
    dq_hom::A
    dqv::A
    dn_in::A
    dn_het::A
    dq_het::A
end

function IceIncrements(namelists::Namelists, domain::Domain)::IceIncrements
    (; ice_setup, nucleation) = namelists.ice
    return IceIncrements(domain, ice_setup, nucleation)
end

function IceIncrements(domain::Domain, ice_setup::NoIce, nucleation::AbstractNucleation)::IceIncrements
    dn_hom = zeros(0, 0, 0)
    dq_hom = zeros(0, 0, 0)
    dqv = zeros(0, 0, 0)
    dn_in = zeros(0, 0, 0)
    dn_het = zeros(0, 0, 0)
    dq_het = zeros(0, 0, 0)

    return IceIncrements(dn_hom, dq_hom, dqv, dn_in, dn_het, dq_het)
end

function IceIncrements(domain::Domain, ice_setup::AbstractIce, nucleation::BothNucleation)::IceIncrements
    (; nxx, nyy, nzz) = domain

    dn_hom = zeros(nxx, nyy, nzz)
    dq_hom = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)
    dn_in = zeros(nxx, nyy, nzz)
    dn_het = zeros(nxx, nyy, nzz)
    dq_het = zeros(nxx, nyy, nzz)

    return IceIncrements(dn_hom, dq_hom, dqv, dn_in, dn_het, dq_het)
end

function IceIncrements(domain::Domain, ice_setup::AbstractIce, nucleation::HomogeneousOnly)::IceIncrements
    (; nxx, nyy, nzz) = domain
    #=
    dn_hom = zeros(nxx, nyy, nzz)
    dq_hom = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)
    dn_in = zeros(nxx, nyy, nzz)
    dn_het = zeros(nxx, nyy, nzz)
    dq_het = zeros(nxx, nyy, nzz)
    =#
    dn_hom = zeros(nxx, nyy, nzz)
    dq_hom = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)
    dn_in = zeros(0, 0, 0)
    dn_het = zeros(0, 0, 0)
    dq_het = zeros(0, 0, 0)
    

    return IceIncrements(dn_hom, dq_hom, dqv, dn_in, dn_het, dq_het)
end

function IceIncrements(domain::Domain, ice_setup::AbstractIce, nucleation::HeterogeneousOnly)::IceIncrements
    (; nxx, nyy, nzz) = domain
    #=
    dn_hom = zeros(nxx, nyy, nzz)
    dq_hom = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)
    dn_in = zeros(nxx, nyy, nzz)
    dn_het = zeros(nxx, nyy, nzz)
    dq_het = zeros(nxx, nyy, nzz)
    =#
    dn_hom = zeros(0, 0, 0)
    dq_hom = zeros(0, 0, 0)
    dqv = zeros(nxx, nyy, nzz)
    dn_in = zeros(nxx, nyy, nzz)
    dn_het = zeros(nxx, nyy, nzz)
    dq_het = zeros(nxx, nyy, nzz)
    

    return IceIncrements(dn_hom, dq_hom, dqv, dn_in, dn_het, dq_het)
end
