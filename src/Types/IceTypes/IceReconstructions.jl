"""
```julia
IceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of ice variables.

The first three dimensions represent physical space, the fourth dimension represents the direction in which the reconstruction was performed and the fifth dimension represents the two cell edges of the reconstruction.

```julia
IceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::IceReconstructions
```

Construct an `IceReconstructions` instance with dimensions depending on the general ice-physics configuration, by dispatching to the appropriate method.

```julia
IceReconstructions(domain::Domain, ice_setup::NoIce)::IceReconstructions
```

Construct an `IceReconstructions` instance with zero-size arrays for configurations without ice physics.

```julia
IceReconstructions(
    domain::Domain,
    ice_setup::AbstractIce,
)::IceReconstructions
```

Construct an `IceReconstructions` instance with zero-initialized arrays.

# Fields

  - `ntilde::A`: Reconstructions of the ice-crystal number concentration.

  - `qtilde::A`: Reconstructions of the ice mixing ratio.

  - `qvtilde::A`: Reconstructions of the water-vapor mixing ratio.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `ice_setup`: General ice-physics configuration.
"""
struct IceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    n_homtilde::A
    q_homtilde::A
    qvtilde::A
    n_intilde::A
    n_hettilde::A
    q_hettilde::A
end

function IceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::IceReconstructions
    (; ice_setup, nucleation) = namelists.ice

    return IceReconstructions(domain, ice_setup, nucleation)
end

function IceReconstructions(domain::Domain, ice_setup::NoIce, nucleation::AbstractNucleation)::IceReconstructions
    n_homtilde = zeros(0, 0, 0, 0, 0)
    q_homtilde = zeros(0, 0, 0, 0, 0)
    qvtilde = zeros(0, 0, 0, 0, 0)
    n_intilde = zeros(0, 0, 0, 0, 0)
    n_hettilde = zeros(0, 0, 0, 0, 0)
    q_hettilde = zeros(0, 0, 0, 0, 0)

    return IceReconstructions(n_homtilde, q_homtilde, qvtilde, n_intilde, n_hettilde, q_hettilde)
end

function IceReconstructions(
    domain::Domain,
    ice_setup::AbstractIce,
    nucleation::BothNucleation,
)::IceReconstructions
    (; nxx, nyy, nzz) = domain

    n_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    q_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    qvtilde = zeros(nxx, nyy, nzz, 3, 2)
    n_intilde = zeros(nxx, nyy, nzz, 3, 2)
    n_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    q_hettilde = zeros(nxx, nyy, nzz, 3, 2)

    return IceReconstructions(n_homtilde, q_homtilde, qvtilde, n_intilde, n_hettilde, q_hettilde)
end

function IceReconstructions(
    domain::Domain,
    ice_setup::AbstractIce,
    nucleation::HeterogeneousOnly,
)::IceReconstructions
    (; nxx, nyy, nzz) = domain

    n_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    q_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    qvtilde = zeros(nxx, nyy, nzz, 3, 2)
    n_intilde = zeros(nxx, nyy, nzz, 3, 2)
    n_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    q_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    #=
    n_homtilde = zeros(0, 0, 0, 0, 0)
    q_homtilde = zeros(0, 0, 0, 0, 0)
    qvtilde = zeros(nxx, nyy, nzz, 3, 2)
    n_intilde = zeros(nxx, nyy, nzz, 3, 2)
    n_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    q_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    =#

    return IceReconstructions(n_homtilde, q_homtilde, qvtilde, n_intilde, n_hettilde, q_hettilde)
end

function IceReconstructions(
    domain::Domain,
    ice_setup::AbstractIce,
    nucleation::HomogeneousOnly,
)::IceReconstructions
    (; nxx, nyy, nzz) = domain

    n_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    q_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    qvtilde = zeros(nxx, nyy, nzz, 3, 2)
    n_intilde = zeros(nxx, nyy, nzz, 3, 2)
    n_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    q_hettilde = zeros(nxx, nyy, nzz, 3, 2)
    #=
    n_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    q_homtilde = zeros(nxx, nyy, nzz, 3, 2)
    qvtilde = zeros(nxx, nyy, nzz, 3, 2)
    n_intilde = zeros(0, 0, 0, 0, 0)
    n_hettilde = zeros(0, 0, 0, 0, 0)
    q_hettilde = zeros(0, 0, 0, 0, 0)
    =#

    return IceReconstructions(n_homtilde, q_homtilde, qvtilde, n_intilde, n_hettilde, q_hettilde)
end
