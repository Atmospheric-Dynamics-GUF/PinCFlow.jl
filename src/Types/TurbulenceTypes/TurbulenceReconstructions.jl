"""
```julia
TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of turbulence energies.

The first three dimensions represent physical space, the fourth represents the physical-space dimension of the reconstruction and the fifth the two directions in which it is computed.

```julia
TurbulenceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with dimensions depending on the general turbulence-transport configuration, by dispatching to the appropriate method.

```julia
TurbulenceReconstructions(
    domain::Domain,
    turbulence_setup::NoTurbulence,
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-size arrays for configurations without turbulence transport.

```julia
TurbulenceReconstructions(
    domain::Domain,
    turbulence_setup::AbstractTurbulence,
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-initialized arrays.

# Fields

  - `chitilde::A`: Reconstructions of a non-dimensional turbulence.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_setup`: General turbulence-transport configuration.
"""
struct TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    chitilde::A
end

function TurbulenceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceReconstructions
    (; turbulence_setup) = namelists.turbulence

    return TurbulenceReconstructions(domain, turbulence_setup)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_setup::NoTurbulence,
)::TurbulenceReconstructions
    return TurbulenceReconstructions(
        [
            zeros(0, 0, 0, 0, 0) for field in fieldnames(TurbulenceReconstructions)
        ]...,
    )
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_setup::AbstractTurbulence,
)::TurbulenceReconstructions
    (; nxx, nyy, nzz) = domain

    return TurbulenceReconstructions(
        [
            zeros(nxx, nyy, nzz, 3, 2) for
            field in fieldnames(TurbulenceReconstructions)
        ]...,
    )
end
