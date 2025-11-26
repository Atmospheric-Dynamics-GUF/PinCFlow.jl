"""
```julia
TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

Arrays for fluxes of turbulence variables.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
TurbulenceFluxes(namelists::Namelists, domain::Domain)::TurbulenceFluxes
```

Construct a `TurbulenceFluxes` instance with dimensions depending on the general turbulence-physics configuration, by dispatching to the appropriate method.

```julia
TurbulenceFluxes(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceFluxes
```

Construct a `TurbulenceFluxes` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulenceFluxes(domain::Domain, turbulence_scheme::TKEScheme)::TurbulenceFluxes
```

Construct a `TurbulenceFluxes` instance with zero-initialized arrays.

# Fields

  - `phitke::A`: Fluxes of the turbulent kinetic energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence-physics configuration.
"""
struct TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phitke::A
end

function TurbulenceFluxes(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceFluxes
    (; turbulence_scheme) = namelists.turbulence

    return TurbulenceFluxes(domain, turbulence_scheme)
end

function TurbulenceFluxes(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceFluxes
    phitke = zeros(0, 0, 0, 0)

    return TurbulenceFluxes(phitke)
end

function TurbulenceFluxes(domain::Domain, turbulence_scheme::TKEScheme)::TurbulenceFluxes
    (; nxx, nyy, nzz) = domain

    phitke = zeros(nxx, nyy, nzz, 3)

    return TurbulenceFluxes(phitke)
end
