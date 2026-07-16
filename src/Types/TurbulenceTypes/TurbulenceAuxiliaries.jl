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
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance with zero-size arrays for configurations without turbulence parameterization.

```julia 
TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance with zero-initialized arrays.

# Fields

  - `shear_production::A`: Contribution of turbulence production due to shear. 

  - `buoyancy_production::A`: Contribution of turbulence production due to the buoyancy.

  - `gw_shear::A`: Gravity-wave induced shear field.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.
"""
struct TurbulenceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    shear_production::A
    buoyancy_production::A
    gw_shear::A
end

function TurbulenceAuxiliaries(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceAuxiliaries
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme return TurbulenceAuxiliaries(
        namelists,
        domain,
        Val(turbulence_scheme),
    )
end

function TurbulenceAuxiliaries(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceAuxiliaries
    print("No Turbulence")
    return TurbulenceAuxiliaries([zeros(0, 0, 0) for i in 1:3]...)
end

function TurbulenceAuxiliaries(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceAuxiliaries
    (; wkb_mode) = namelists.wkb
    (; nxx, nyy, nzz) = domain

    shear_production = zeros(nxx, nyy, nzz)
    buoyancy_production = zeros(nxx, nyy, nzz)

    if wkb_mode != :MultiColumn
        gw_shear = zeros(0, 0, 0)
    else
        gw_shear = zeros(nxx, nyy, nzz)
    end

    return TurbulenceAuxiliaries(
        shear_production,
        buoyancy_production,
        gw_shear,
    )
end
