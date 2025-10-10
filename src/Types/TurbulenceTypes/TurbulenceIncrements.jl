"""
```julia
TurbulenceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of turbulence variables.

```julia
TurbulenceIncrements(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with dimensions depending on the general turbulence-physics configuration, by dispatching to the appropriate method.

```julia
TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::AbstractTurbulence,
)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with zero-initialized arrays.

# Fields

  - `dtke::A`: Runge-Kutta update of the turbulent kinetic energy.

  - `dtte::A`: Runge-Kutta update of the total turbulent energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence-physics configuration.
"""
struct TurbulenceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
    dtke::A
    dtte::A
end

function TurbulenceIncrements(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceIncrements
    (; turbulence_scheme) = namelists.turbulence

    return TurbulenceIncrements(domain, turbulence_scheme)
end

function TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceIncrements
    dtke = zeros(0, 0, 0)
    dtte = zeros(0, 0, 0)

    return TurbulenceIncrements(dtke, dtte)
end

function TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::AbstractTurbulence,
)::TurbulenceIncrements
    (; nxx, nyy, nzz) = domain

    dtke = zeros(nxx, nyy, nzz)
    dtte = zeros(nxx, nyy, nzz)

    return TurbulenceIncrements(dtke, dtte)
end
