"""
```julia
TurbulenceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of turbulence variables.

```julia
TurbulenceTendencies(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceTendencies
```

Construct a `TurbulenceTendencies` instance with dimensions depending on the general turbulence-physics configuration, by dispatching to the appropriate method.

```julia
TurbulenceTendencies(
    domain::Domain,
    turbulencesetup::NoTurbulence,
)::TurbulenceTendencies
```

Construct a `TurbulenceTendencies` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulenceTendencies(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)::TurbulenceTendencies
```

Construct a `TurbulenceTendencies` instance with zero-initialized arrays.

# Fields

  - `dtke::A`: Runge-Kutta update of the turbulent kinetic energy.

  - `dtte::A`: Runge-Kutta update of the total turbulent energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulencesetup`: General turbulence-physics configuration.
"""
struct TurbulenceTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dtke::A
    dtte::A
end

function TurbulenceTendencies(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceTendencies
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceTendencies(domain, turbulencesetup)
end

function TurbulenceTendencies(
    domain::Domain,
    turbulencesetup::NoTurbulence,
)::TurbulenceTendencies
    dtke = zeros(0, 0, 0)
    dtte = zeros(0, 0, 0)

    return TurbulenceTendencies(dtke, dtte)
end

function TurbulenceTendencies(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)::TurbulenceTendencies
    (; nxx, nyy, nzz) = domain

    dtke = zeros(nxx, nyy, nzz)
    dtte = zeros(nxx, nyy, nzz)

    return TurbulenceTendencies(dtke, dtte)
end
