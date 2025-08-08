"""
```julia
TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

Arrays for fluxes of turbulence variables.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
TurbulenceFluxes(namelists::Namelists, domain::Domain)
```

Construct a `TurbulenceFluxes` instance with dimensions depending on the general turbulence-physics configuration, by dispatching to the appropriate method.

```julia
TurbulenceFluxes(domain::Domain, turbulencesetup::NoTurbulence)
```

Construct a `TurbulenceFluxes` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulenceFluxes(domain::Domain, turbulencesetup::AbstractTurbulence)
```

Construct a `TurbulenceFluxes` instance with zero-initialized arrays.

# Fields

- `phitke::A`: Fluxes of the turbulent kinetic energy.
- `phitte::A`: Fluxes of the total turbulent energy.

# Arguments

- `namelists`: Namelists with all model parameters.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `turbulencesetup`: General turbulence-physics configuration.
"""
struct TurbulenceFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phitke::A
    phitte::A
end

function TurbulenceFluxes(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceFluxes(domain, turbulencesetup)
end

function TurbulenceFluxes(domain::Domain, turbulencesetup::NoTurbulence)
    phitke = zeros(0, 0, 0, 0)
    phitte = zeros(0, 0, 0, 0)

    return TurbulenceFluxes(phitke, phitte)
end

function TurbulenceFluxes(domain::Domain, turbulencesetup::AbstractTurbulence)
    (; nxx, nyy, nzz) = domain

    phitke = zeros(nxx, nyy, nzz, 3)
    phitte = zeros(nxx, nyy, nzz, 3)

    return TurbulenceFluxes(phitke, phitte)
end
