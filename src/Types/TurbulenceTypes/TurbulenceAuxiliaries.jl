"""
```julia
TurbulenceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Shear and buoyancy production terms.

```julia
TurbulenceAuxiliaries(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance with dimensions depending on the general turbulence parameterization configuration, by dispatching to the appropriate method. 

```julia 
TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance with zero-size arrays for configurations without turbulence parameterization.

```julia 
TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::TKEScheme,
)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance with zero-initialized arrays.

# Fields

  - `shearproduction::A`: Contribution of turbulence production due to shear. 

  - `buoyancyproduction::A`: Contribution of turbulence production due to the buoyancy.

# Arguments

  - `namelists`: Namelists with all model paramters.

  - `domain`: Collection of domain-decomposition and MPI-communication paramters.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
struct TurbulenceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    shearproduction::A
    buoyancyproduction::A
end

function TurbulenceAuxiliaries(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceAuxiliaries
    (; turbulence_scheme) = namelists.turbulence

    return TurbulenceAuxiliaries(domain, turbulence_scheme)
end

function TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceAuxiliaries
    return TurbulenceAuxiliaries([zeros(0, 0, 0) for i in 1:2]...)
end

function TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::TKEScheme,
)::TurbulenceAuxiliaries
    (; nxx, nyy, nzz) = domain
    return TurbulenceAuxiliaries([zeros(nxx, nyy, nzz) for i in 1:2]...)
end
