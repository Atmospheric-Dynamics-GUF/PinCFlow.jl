"""
```julia
TurbulenceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of turbulence variables.

```julia
TurbulenceIncrements(namelists::Namelists, domain::Domain)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with dimensions depending on the turbulence parameterization configuration, by dispatching to the appropriate method.

```julia
TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with zero-size arrays for configurations without turbulence parameterization.

```julia
TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with zero-initialized arrays.

# Fields

  - `dtke::A`: Runge-Kutta update of the turbulent kinetic energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
struct TurbulenceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
    dtke::A
end

function TurbulenceIncrements(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceIncrements
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme return TurbulenceIncrements(domain, Val(turbulence_scheme))
end

function TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceIncrements
    dtke = zeros(0, 0, 0)

    return TurbulenceIncrements(dtke)
end

function TurbulenceIncrements(
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceIncrements
    (; nxx, nyy, nzz) = domain

    dtke = zeros(nxx, nyy, nzz)

    return TurbulenceIncrements(dtke)
end
