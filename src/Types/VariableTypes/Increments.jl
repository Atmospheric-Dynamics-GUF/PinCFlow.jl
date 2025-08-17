"""
```julia
Increments{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Container for the Runge-Kutta updates of the prognostic variables, as well as the Exner-pressure update of the Poisson solver.

```julia
Increments(namelists::Namelists, domain::Domain)::Increments
```

Create a `Increments` instance with dimensions depending on whether or not the model is compressible, by dispatching to the appropriate method.

```julia
Increments(domain::Domain, model::AbstractModel)::Increments
```

Create a `Increments` instance in non-compressible modes, with a zero-size array for the mass-weighted potential-temperature update.

```julia
Increments(domain::Domain, model::Compressible)::Increments
```

Create a `Increments` instance in compressible mode.

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
struct Increments{
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

function Increments(namelists::Namelists, domain::Domain)::Increments
    (; model) = namelists.setting
    return Increments(domain, model)
end

function Increments(domain::Domain, model::AbstractModel)::Increments
    (; nxx, nyy, nzz) = domain

    # Initialize the increments.
    (drho, drhop, du, dv, dw, dpip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    dp = zeros(0, 0, 0)

    # Return a Variables instance.
    return Increments(drho, drhop, du, dv, dw, dpip, dp)
end

function Increments(domain::Domain, model::Compressible)::Increments
    (; nxx, nyy, nzz) = domain

    # Initialize the increments.
    (drho, drhop, du, dv, dw, dpip, dp) = (zeros(nxx, nyy, nzz) for i in 1:7)

    # Return a Variables instance.
    return Increments(drho, drhop, du, dv, dw, dpip, dp)
end
