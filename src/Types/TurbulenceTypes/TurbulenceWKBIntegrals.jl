"""
```julia
TurbulenceWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Gravity-wave shear field.

```julia 
TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance by dispatching to the appropriate method.

```julia
TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance with zero-size arrays for configurations without turbulence parameterization.

```julia 
TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance by dispatching to the appropriate method.

```julia 
TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance with zero-size arrays for non-WKB configurations.

```julia
TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance with zero-initialized arrays.

# Fields 

  - `gw_shear::A`: Gravity-wave shear field.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TurbulenceWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    gw_shear::A
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBIntegrals
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme return TurbulenceWKBIntegrals(
        namelists,
        domain,
        Val(turbulence_scheme),
    )
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceWKBIntegrals
    return TurbulenceWKBIntegrals([zeros(0, 0, 0) for i in 1:1]...)
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceWKBIntegrals
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TurbulenceWKBIntegrals(
        namelists,
        domain,
        Val(wkb_mode),
    )
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TurbulenceWKBIntegrals
    return TurbulenceWKBIntegrals([zeros(0, 0, 0) for i in 1:1]...)
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBIntegrals
    (; nxx, nyy, nzz) = domain

    return TurbulenceWKBIntegrals([zeros(nxx, nyy, nzz) for i in 1:1]...)
end
