"""
```julia
Fluxes{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
}
```

Arrays for fluxes needed in the computation of the left-hand sides.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
Fluxes(namelists::Namelists, domain::Domain)::Fluxes
```

Construct a `Fluxes` instance with dimensions depending on whether or not the model is compressible, by dispatching to the appropriate method.

```julia
Fluxes(domain::Domain, model::Union{Boussinesq, PseudoIncompressible})::Fluxes
```

Construct a `Fluxes` instance in non-compressible modes, with a zero-size array for mass-weighted potential-temperature fluxes.

```julia
Fluxes(domain::Domain, model::Compressible)::Fluxes
```

Construct a `Fluxes` instance in compressible mode.

# Fields

  - `phirho::A`: Density fluxes.

  - `phirhop::A`: Density-fluctuations fluxes.

  - `phiu::A`: Zonal-momentum fluxes.

  - `phiv::A`: Meridional-momentum fluxes.

  - `phiw::A`: Transformed-vertical-momentum fluxes.

  - `phip::B`: Mass-weighted potential-temperature fluxes.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `model`: Dynamic equations.
"""
struct Fluxes{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
}
    phirho::A
    phirhop::A
    phiu::A
    phiv::A
    phiw::A
    phitheta::A
    phip::B
end

function Fluxes(namelists::Namelists, domain::Domain)::Fluxes
    (; model) = namelists.atmosphere
    return Fluxes(domain, model)
end

function Fluxes(
    domain::Domain,
    model::Union{Boussinesq, PseudoIncompressible},
)::Fluxes
    (; nxx, nyy, nzz) = domain

    return Fluxes([zeros(nxx, nyy, nzz, 3) for i in 1:6]..., zeros(0, 0, 0, 0))
end

function Fluxes(domain::Domain, model::Compressible)::Fluxes
    (; nxx, nyy, nzz) = domain

    return Fluxes([zeros(nxx, nyy, nzz, 3) for i in 1:7]...)
end
