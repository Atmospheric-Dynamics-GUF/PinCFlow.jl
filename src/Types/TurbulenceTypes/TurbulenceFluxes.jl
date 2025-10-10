"""
```julia
TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

Arrays for fluxes of turbulence energies.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
TurbulenceFluxes(namelists::Namelists, domain::Domain)::TurbulenceFluxes
```

Construct a `TurbulenceFluxes` instance with dimensions depending on the general turbulence-transport configuration, by dispatching to the appropriate method.

```julia
TurbulenceFluxes(domain::Domain, turbulence_setup::NoTurbulence)::TurbulenceFluxes
```

Construct a `TurbulenceFluxes` instance with zero-size arrays for configurations without turbulence transport.

```julia
TurbulenceFluxes(domain::Domain, turbulence_setup::AbstractTurbulence)::TurbulenceFluxes
```

Construct a `TurbulenceFluxes` instance with zero-initialized arrays.

# Fields

  - `phichi::A`: Fluxes of a non-dimensional turbulence.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_setup`: General turbulence-transport configuration.
"""
struct TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phichi::A
end

function TurbulenceFluxes(namelists::Namelists, domain::Domain)::TurbulenceFluxes
    (; turbulence_setup) = namelists.turbulence

    return TurbulenceFluxes(domain, turbulence_setup)
end

function TurbulenceFluxes(domain::Domain, turbulence_setup::NoTurbulence)::TurbulenceFluxes
    return TurbulenceFluxes(
        [zeros(0, 0, 0, 0) for field in fieldnames(TurbulenceFluxes)]...,
    )
end

function TurbulenceFluxes(
    domain::Domain,
    turbulence_setup::AbstractTurbulence,
)::TurbulenceFluxes
    (; nxx, nyy, nzz) = domain

    return TurbulenceFluxes(
        [zeros(nxx, nyy, nzz, 3) for field in fieldnames(TurbulenceFluxes)]...,
    )
end
