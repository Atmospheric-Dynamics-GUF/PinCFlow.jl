"""
```julia
TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of turbulence variables.

The first three dimensions represent physical space, the fourth dimension represents the direction in which the reconstruction was performed and the fifth dimension represents the two cell edges of the reconstruction.

```julia
TurbulenceReconstructions(namelists::Namelists, domain::Domain)
```

Construct a `TurbulenceReconstructions` instance with dimensions depending on the general turbulence-physics configuration, by dispatching to the appropriate method.

```julia
TurbulenceReconstructions(domain::Domain, turbulencesetup::NoTurbulence)
```

Construct a `TurbulenceReconstructions` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulenceReconstructions(domain::Domain, turbulencesetup::AbstractTurbulence)
```

Construct a `TurbulenceReconstructions` instance with zero-initialized arrays.

# Fields

- `tketilde::A`: Reconstructions of the turbulent kinetic energy.
- `ttetilde::A`: Reconstructions of the total turbulent energy.

# Arguments

- `namelists`: Namelists with all model parameters.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `turbulencesetup`: General turbulence-physics configuration.
"""
struct TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    tketilde::A
    ttetilde::A
end

function TurbulenceReconstructions(namelists::Namelists, domain::Domain)
    (; turbulencesetup) = namelists.turbulence

    return TurbulenceReconstructions(domain, turbulencesetup)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulencesetup::NoTurbulence,
)
    tketilde = zeros(0, 0, 0, 0, 0)
    ttetilde = zeros(0, 0, 0, 0, 0)

    return TurbulenceReconstructions(tketilde, ttetilde)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulencesetup::AbstractTurbulence,
)
    (; nxx, nyy, nzz) = domain

    tketilde = zeros(nxx, nyy, nzz, 3, 2)
    ttetilde = zeros(nxx, nyy, nzz, 3, 2)

    return TurbulenceReconstructions(tketilde, ttetilde)
end
