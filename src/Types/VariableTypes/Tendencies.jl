"""
```julia
Tendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Container for the Runge-Kutta updates of the prognostic variables, as well as the Exner-pressure update of the Poisson solver.

```julia
Tendencies(namelists::Namelists, domain::Domain)
```

Create a `Tendencies` instance with dimensions depending on whether or not the model is compressible, by dispatching to the appropriate method.

```julia
Tendencies(domain::Domain, model::AbstractModel)
```

Create a `Tendencies` instance in non-compressible modes, with a zero-size array for the mass-weighted potential-temperature update.

```julia
Tendencies(domain::Domain, model::Compressible)
```

Create a `Tendencies` instance in compressible mode.

# Fields

  - `drho::A`: Density update.

  - `drhop::A`: Density-fluctuation update.

  - `du::A`: Zonal-momentum update.

  - `dv::A`: Meridional-momentum update.

  - `dw::A`: Transformed-vertical-momentum update.

  - `dpip::A`: Exner-pressure update.

  - `dp::B`: Mass-weighted potential-temperature update.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `model`: Dynamic equations.
"""
struct Tendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    drho::A
    drhop::A
    du::A
    dv::A
    dw::A
    dpip::A
    dp::B
end

function Tendencies(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Tendencies(domain, model)
end

function Tendencies(domain::Domain, model::AbstractModel)
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    dp = zeros(0, 0, 0)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip, dp)
end

function Tendencies(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip, dp) = (zeros(nxx, nyy, nzz) for i in 1:7)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip, dp)
end
