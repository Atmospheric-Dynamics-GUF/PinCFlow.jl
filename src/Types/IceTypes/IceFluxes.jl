"""
```julia
IceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

Arrays for fluxes of ice variables.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
IceFluxes(namelists::Namelists, domain::Domain)::IceFluxes
```

Construct an `IceFluxes` instance with dimensions depending on the general ice-physics configuration, by dispatching to the appropriate method.

```julia
IceFluxes(domain::Domain, ice_setup::NoIce)::IceFluxes
```

Construct an `IceFluxes` instance with zero-size arrays for configurations without ice physics.

```julia
IceFluxes(domain::Domain, ice_setup::AbstractIce)::IceFluxes
```

Construct an `IceFluxes` instance with zero-initialized arrays.

# Fields

  - `phin::A`: Fluxes of the ice-crystal number concentration.

  - `phiq::A`: Fluxes of the ice mixing ratio.

  - `phiqv::A`: Fluxes of the water-vapor mixing ratio.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `ice_setup`: General ice-physics configuration.
"""
struct IceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phin_hom::A
    phiq_hom::A
    phiqv::A
    phin_in::A
    phin_het::A
    phiq_het::A
end

function IceFluxes(namelists::Namelists, domain::Domain)::IceFluxes
    (; ice_setup, nucleation) = namelists.ice

    return IceFluxes(domain, ice_setup, nucleation)
end

function IceFluxes(domain::Domain, ice_setup::NoIce, nucleation::AbstractNucleation)::IceFluxes
    phin_hom = zeros(0, 0, 0, 0)
    phiq_hom = zeros(0, 0, 0, 0)
    phiqv = zeros(0, 0, 0, 0)
    phin_in = zeros(0, 0, 0, 0)
    phin_het = zeros(0, 0, 0 ,0)
    phiq_het = zeros(0, 0, 0, 0)

    return IceFluxes(phin_hom, phiq_hom, phiqv, phin_in, phin_het, phiq_het)
end

# homogeneous nucleation
function IceFluxes(domain::Domain, ice_setup::AbstractIce, nucleation::HomogeneousOnly)::IceFluxes
    (; nxx, nyy, nzz) = domain
    #=
    phin_hom = zeros(nxx, nyy, nzz, 3)
    phiq_hom = zeros(nxx, nyy, nzz, 3)
    phiqv = zeros(nxx, nyy, nzz, 3)
    phin_in = zeros(nxx, nyy, nzz, 3)
    phin_het = zeros(nxx, nyy, nzz, 3)
    phiq_het = zeros(nxx, nyy, nzz, 3)
    =#
    phin_hom = zeros(nxx, nyy, nzz, 3)
    phiq_hom = zeros(nxx, nyy, nzz, 3)
    phiqv = zeros(nxx, nyy, nzz, 3)
    phin_in = zeros(0, 0, 0, 0)
    phin_het = zeros(0, 0, 0, 0)
    phiq_het = zeros(0, 0, 0, 0)
    

    return IceFluxes(phin_hom, phiq_hom, phiqv, phin_in, phin_het, phiq_het)
end

# heterogeneous nucleation
function IceFluxes(domain::Domain, ice_setup::AbstractIce, nucleation::HeterogeneousOnly)::IceFluxes
    (; nxx, nyy, nzz) = domain
    #=
    phin_hom = zeros(nxx, nyy, nzz, 3)
    phiq_hom = zeros(nxx, nyy, nzz, 3)
    phiqv = zeros(nxx, nyy, nzz, 3)
    phin_in = zeros(nxx, nyy, nzz, 3)
    phin_het = zeros(nxx, nyy, nzz, 3)
    phiq_het = zeros(nxx, nyy, nzz, 3)
    =#
    phin_hom = zeros(0, 0, 0, 0)
    phiq_hom = zeros(0, 0, 0, 0)
    phiqv = zeros(nxx, nyy, nzz, 3)
    phin_in = zeros(nxx, nyy, nzz, 3)
    phin_het = zeros(nxx, nyy, nzz, 3)
    phiq_het = zeros(nxx, nyy, nzz, 3)
    

    return IceFluxes(phin_hom, phiq_hom, phiqv, phin_in, phin_het, phiq_het)
end

# both nucleation
function IceFluxes(domain::Domain, ice_setup::AbstractIce, nucleation::BothNucleation)::IceFluxes
    (; nxx, nyy, nzz) = domain

    phin_hom = zeros(nxx, nyy, nzz, 3)
    phiq_hom = zeros(nxx, nyy, nzz, 3)
    phiqv = zeros(nxx, nyy, nzz, 3)
    phin_in = zeros(nxx, nyy, nzz, 3)
    phin_het = zeros(nxx, nyy, nzz, 3)
    phiq_het = zeros(nxx, nyy, nzz, 3)

    return IceFluxes(phin_hom, phiq_hom, phiqv, phin_in, phin_het, phiq_het)
end
