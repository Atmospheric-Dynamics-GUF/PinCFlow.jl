"""
```julia
Reconstructions{
    A <: AbstractArray{<:AbstractFloat, 5},
    B <: AbstractArray{<:AbstractFloat, 5},
    C <: AbstractArray{<:AbstractFloat, 5},
    D <: AbstractArray{<:AbstractFloat, 5},
    E <: AbstractArray{<:AbstractFloat, 5},
}
```

Arrays for the reconstructions of prognostic variables.

The first three dimensions represent physical space, the fourth represents the physical-space dimension of the reconstruction and the fifth the two directions in which it is computed.

```julia
Reconstructions(namelists::Namelists, domain::Domain)::Reconstructions
```

Construct a `Reconstructions` instance with dimensions depending on whether or not the model is Boussinesq, by dispatching to the appropriate method.

```julia
Reconstructions(domain::Domain, model::Boussinesq)::Reconstructions
```

Construct a `Reconstructions` instance in Boussinesq mode, with a zero-size array for density reconstructions.

```julia
Reconstructions(
    domain::Domain,
    model::Union{PseudoIncompressible, Compressible},
)::Reconstructions
```

Construct a `Reconstructions` instance in non-Boussinesq modes.

# Fields

  - `rhotilde::A`: Reconstructed density.

  - `rhoptilde::B`: Reconstructed density fluctuations.

  - `utilde::C`: Reconstructed zonal momentum.

  - `vtilde::D`: Reconstructed meridional momentum.

  - `wtilde::E`: Reconstructed vertical momentum.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `model`: Dynamic equations.
"""
struct Reconstructions{
    A <: AbstractArray{<:AbstractFloat, 5},
    B <: AbstractArray{<:AbstractFloat, 5},
    C <: AbstractArray{<:AbstractFloat, 5},
    D <: AbstractArray{<:AbstractFloat, 5},
    E <: AbstractArray{<:AbstractFloat, 5},
}
    rhotilde::A
    rhoptilde::B
    utilde::C
    vtilde::D
    wtilde::E
end

function Reconstructions(namelists::Namelists, domain::Domain)::Reconstructions
    (; model) = namelists.atmosphere

    return Reconstructions(domain, model)
end

function Reconstructions(domain::Domain, model::Boussinesq)::Reconstructions
    (; nxx, nyy, nzz) = domain

    return Reconstructions(
        zeros(0, 0, 0, 0, 0),
        [zeros(nxx, nyy, nzz, 3, 2) for i in 1:4]...,
    )
end

function Reconstructions(
    domain::Domain,
    model::Union{PseudoIncompressible, Compressible},
)::Reconstructions
    (; nxx, nyy, nzz) = domain

    return Reconstructions([zeros(nxx, nyy, nzz, 3, 2) for i in 1:5]...)
end
