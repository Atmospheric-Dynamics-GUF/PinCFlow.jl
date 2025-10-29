"""
```julia
Fluxes{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
    C <: AbstractArray{<:AbstractFloat, 4},
    D <: AbstractArray{<:AbstractFloat, 4},
    E <: AbstractArray{<:AbstractFloat, 4},
    F <: AbstractArray{<:AbstractFloat, 4},
    G <: AbstractArray{<:AbstractFloat, 4},
}
```

Arrays for fluxes needed in the computation of the left-hand sides.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
Fluxes(namelists::Namelists, domain::Domain)::Fluxes
```

Construct a `Fluxes` instance with dimensions depending on whether or not the model is compressible, by dispatching to the appropriate method.

```julia
Fluxes(domain::Domain, model::Boussinesq)::Fluxes
```

Construct a `Fluxes` instance in Boussinesq mode, with zero-size arrays for the density and mass-weighted potential-temperature fluxes.

```julia
Fluxes(domain::Domain, model::PseudoIncompressible)::Fluxes
```

Construct a `Fluxes` instance in pseudo-incompressible mode, with a zero-size array for mass-weighted potential-temperature fluxes.

```julia
Fluxes(domain::Domain, model::Compressible)::Fluxes
```

Construct a `Fluxes` instance in compressible mode.

# Fields

  - `phirho::A`: Density fluxes.

  - `phirhop::B`: Density-fluctuations fluxes.

  - `phiu::C`: Zonal-momentum fluxes.

  - `phiv::D`: Meridional-momentum fluxes.

  - `phiw::E`: Transformed-vertical-momentum fluxes.

  - `phitheta::F`: Potential temperature fluxes.

  - `phip::G`: Mass-weighted potential-temperature fluxes.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `model`: Dynamic equations.
"""
struct Fluxes{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
    C <: AbstractArray{<:AbstractFloat, 4},
    D <: AbstractArray{<:AbstractFloat, 4},
    E <: AbstractArray{<:AbstractFloat, 4},
    F <: AbstractArray{<:AbstractFloat, 4},
    G <: AbstractArray{<:AbstractFloat, 4},
}
    phirho::A
    phirhop::B
    phiu::C
    phiv::D
    phiw::E
    phitheta::F
    phip::G
end

function Fluxes(namelists::Namelists, domain::Domain)::Fluxes
    (; model) = namelists.atmosphere
    return Fluxes(domain, model)
end

function Fluxes(domain::Domain, model::Boussinesq)::Fluxes
    (; nxx, nyy, nzz) = domain

    return Fluxes(
        zeros(0, 0, 0, 0),
        [zeros(nxx, nyy, nzz, 3) for i in 1:5]...,
        zeros(0, 0, 0, 0),
    )
end

function Fluxes(domain::Domain, model::PseudoIncompressible)::Fluxes
    (; nxx, nyy, nzz) = domain

    return Fluxes([zeros(nxx, nyy, nzz, 3) for i in 1:6]..., zeros(0, 0, 0, 0))
end

function Fluxes(domain::Domain, model::Compressible)::Fluxes
    (; nxx, nyy, nzz) = domain

    return Fluxes([zeros(nxx, nyy, nzz, 3) for i in 1:7]...)
end
