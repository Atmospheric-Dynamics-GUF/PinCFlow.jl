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
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-size arrays for configurations without turbulence parameterization.

```julia
TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-initialized arrays.

```julia
TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::Val{:TTEScheme},
)::TurbulenceReconstructions
```

Construct a `TurbulenceReconstructions` instance with zero-initialized arrays.

# Fields

  - `tketilde::A`: Reconstructions of the non-dimensional turbulent kinetic energy.

  - `tapetilde::A`: Reconstructions of the non-dimensional turbulent available potential energy.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
struct TurbulenceReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    tketilde::A
    tapetilde::A
end

function TurbulenceReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceReconstructions
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme TurbulenceReconstructions(
        domain,
        Val(turbulence_scheme),
    )
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceReconstructions
    TurbulenceReconstructions((zeros(0, 0, 0, 0, 0) for i in 1:2)...)
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceReconstructions
    (; nxx, nyy, nzz) = domain

    TurbulenceReconstructions(zeros(nxx, nyy, nzz, 3, 2), zeros(0, 0, 0, 0, 0))
end

function TurbulenceReconstructions(
    domain::Domain,
    turbulence_scheme::Val{:TTEScheme},
)::TurbulenceReconstructions
    (; nxx, nyy, nzz) = domain

    TurbulenceReconstructions((zeros(nxx, nyy, nzz, 3, 2) for i in 1:2)...)
end
