"""
```julia
TurbulenceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of turbulence energies.

```julia
TurbulenceIncrements(namelists::Namelists, domain::Domain)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with dimensions depending on the general turbulence-transport configuration, by dispatching to the appropriate method.

```julia
TurbulenceIncrements(domain::Domain, turbulence_setup::NoTurbulence)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with zero-size arrays for configurations without turbulence transport.

```julia
TurbulenceIncrements(domain::Domain, turbulence_setup::AbstractTurbulence)::TurbulenceIncrements
```

Construct a `TurbulenceIncrements` instance with zero-initialized arrays.

# Fields

  - `dchi::A`: Runge-Kutta update of a non-dimensional turbulence.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_setup`: General turbulence-transport configuration.
"""
struct TurbulenceIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
    dchi::A
end

function TurbulenceIncrements(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceIncrements
    (; turbulence_setup) = namelists.turbulence
    return TurbulenceIncrements(domain, turbulence_setup)
end

function TurbulenceIncrements(
    domain::Domain,
    turbulence_setup::NoTurbulence,
)::TurbulenceIncrements
    return TurbulenceIncrements(
        [zeros(0, 0, 0) for field in fieldnames(TurbulenceIncrements)]...,
    )
end

function TurbulenceIncrements(
    domain::Domain,
    turbulence_setup::AbstractTurbulence,
)::TurbulenceIncrements
    (; nxx, nyy, nzz) = domain

    return TurbulenceIncrements(
        [zeros(nxx, nyy, nzz) for field in fieldnames(TurbulenceIncrements)]...,
    )
end
