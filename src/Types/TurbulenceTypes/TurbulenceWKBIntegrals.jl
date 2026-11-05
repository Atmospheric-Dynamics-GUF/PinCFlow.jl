"""
```julia
TurbulenceWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Integrals of gravity-wave quantities that impact turbulence.

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
    turbulence_scheme::Union{
        Val{:NoTurbulence},
        Val{:TKEScheme},
        Val{:TTEScheme},
    },
    wkb_mode::Union{
        Val{:NoWKB},
        Val{:SteadyState},
        Val{:SingleColumn},
        Val{:MultiColumn},
    },
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance with zero-size arrays for non-WKB or no turbulence parameterization configurations.

```julia
TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Union{Val{:TKEScheme}, Val{:TTEScheme}},
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBIntegrals
```

Construct a `TurbulenceWKBIntegrals` instance with zero-initialized arrays if `state.namelists.turbulence.gw_coupling == true`, otherwise the arrays are zero-size.

# Fields

  - `gw_shear::A`: Turbulence production due to gravity-wave shear.

  - `gw_convective::A`: Turbulence production due to gravity-wave convection.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TurbulenceWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    gw_shear::A
    gw_convective::A
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBIntegrals
    (; turbulence_scheme) = namelists.turbulence
    (; wkb_mode) = namelists.wkb

    @dispatch_turbulence_scheme @dispatch_wkb_mode TurbulenceWKBIntegrals(
        namelists,
        domain,
        Val(turbulence_scheme),
        Val(wkb_mode),
    )
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Union{
        Val{:NoTurbulence},
        Val{:TKEScheme},
        Val{:TTEScheme},
    },
    wkb_mode::Union{
        Val{:NoWKB},
        Val{:SteadyState},
        Val{:SingleColumn},
        Val{:MultiColumn},
    },
)::TurbulenceWKBIntegrals
    TurbulenceWKBIntegrals((zeros(0, 0, 0) for i in 1:2)...)
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBIntegrals
    (; nxx, nyy, nzz) = domain
    (; gw_coupling) = namelists.turbulence

    if gw_coupling
        TurbulenceWKBIntegrals(zeros(nxx, nyy, nzz), zeros(0, 0, 0))
    else
        TurbulenceWKBIntegrals((zeros(0, 0, 0) for i in 1:2)...)
    end
end

function TurbulenceWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TTEScheme},
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBIntegrals
    (; nxx, nyy, nzz) = domain
    (; gw_coupling) = namelists.turbulence

    if gw_coupling
        TurbulenceWKBIntegrals((zeros(nxx, nyy, nzz) for i in 1:2)...)
    else
        TurbulenceWKBIntegrals((zeros(0, 0, 0) for i in 1:2)...)
    end
end
