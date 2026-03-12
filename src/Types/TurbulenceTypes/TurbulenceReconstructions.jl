"""
```julia
TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of turbulence variables.

The first three dimensions represent physical space, the fourth represents the physical-space dimension of the reconstruction and the fifth the two directions in which it is computed.

```julia
TurbulenceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with dimensions depending on the general turbulence parameterization configuration, by dispatching to the appropriate method.

```julia
TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-size arrays for configurations without turbulence parameterization.

```julia
TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::TKEScheme,
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-initialized arrays.

# Fields

  - `tketilde::A`: Reconstructions of the non-dimensional turbulent kinetic energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
struct TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    tketilde::A
end

function TurbulenceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceReconstructions
    (; turbulence_scheme) = namelists.turbulence

    return TurbulenceReconstructions(domain, turbulence_scheme)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceReconstructions
    tketilde = zeros(0, 0, 0, 0, 0)

    return TurbulenceReconstructions(tketilde)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::TKEScheme,
)::TurbulenceReconstructions
    (; nxx, nyy, nzz) = domain

    tketilde = zeros(nxx, nyy, nzz, 3, 2)

    return TurbulenceReconstructions(tketilde)
end
