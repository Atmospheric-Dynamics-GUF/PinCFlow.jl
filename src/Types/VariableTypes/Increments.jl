"""
```julia
Increments{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
```

Container for the Runge-Kutta updates of the prognostic variables, as well as the Exner-pressure update of the Poisson solver.

```julia
Increments(namelists::Namelists, domain::Domain)::Increments
```

Create an `Increments` instance with dimensions depending on the model configuration, by dispatching to the appropriate method.

```julia
Increments(domain::Domain, model::Boussinesq)::Increments
```

Create an `Increments` instance in Boussinesq mode, with zero-size arrays for the density and mass-weighted potential-temperature update.

```julia
Increments(domain::Domain, model::PseudoIncompressible)::Increments
```

Create an `Increments` instance in pseudo-incompressible mode, with a zero-size array for the mass-weighted potential-temperature update.

```julia
Increments(domain::Domain, model::Compressible)::Increments
```

Create an `Increments` instance in compressible mode.

# Fields

  - `drho::A`: Density update.

  - `drhop::B`: Density-fluctuations update.

  - `du::B`: Zonal-momentum update.

  - `dv::B`: Meridional-momentum update.

  - `dw::B`: Transformed-vertical-momentum update.

  - `dpip::B`: Exner-pressure update.

  - `dp::C`: Mass-weighted potential-temperature update.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `model`: Dynamic equations.
"""
struct Increments{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
    drho::A
    drhop::B
    du::B
    dv::B
    dw::B
    dpip::B
    dp::C
end

function Increments(namelists::Namelists, domain::Domain)::Increments
    (; model) = namelists.atmosphere

    return Increments(domain, model)
end

function Increments(domain::Domain, model::Boussinesq)::Increments
    (; nxx, nyy, nzz) = domain

    return Increments(
        zeros(0, 0, 0),
        [zeros(nxx, nyy, nzz) for i in 1:5]...,
        zeros(0, 0, 0),
    )
end

function Increments(domain::Domain, model::PseudoIncompressible)::Increments
    (; nxx, nyy, nzz) = domain

    return Increments([zeros(nxx, nyy, nzz) for i in 1:6]..., zeros(0, 0, 0))
end

function Increments(domain::Domain, model::Compressible)::Increments
    (; nxx, nyy, nzz) = domain

    return Increments([zeros(nxx, nyy, nzz) for i in 1:7]...)
end
