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
Fluxes(namelists::Namelists, domain::Domain)
```

Construct a `Fluxes` instance with dimensions depending on whether or not the model is compressible, by dispatching to the appropriate method.

```julia
Fluxes(domain::Domain, model::AbstractModel)
```

Construct a `Fluxes` instance in non-compressible modes, with a zero-size array for mass-weighted potential-temperature fluxes.

```julia
Fluxes(domain::Domain, model::Compressible)
```

Construct a `Fluxes` instance in compressible mode.

# Fields

- `phirho::A`: Density fluxes.
- `phirhop::A`: Density-fluctuation fluxes.
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
    phip::B
end

function Fluxes(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Fluxes(domain, model)
end

function Fluxes(domain::Domain, model::AbstractModel)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw) = (zeros(nxx, nyy, nzz, 3) for i in 1:5)
    phip = zeros(0, 0, 0, 0)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end

function Fluxes(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw, phip) =
        (zeros(nxx, nyy, nzz, 3) for i in 1:6)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end
