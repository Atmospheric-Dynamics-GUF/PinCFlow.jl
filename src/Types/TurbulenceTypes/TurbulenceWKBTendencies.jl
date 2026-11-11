"""
```julia
TurbulenceWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Turbulence impact of unresolved gravity waves.

```julia
TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance by dispatching to the appropriate method.

```julia
TurbulenceWKBTendencies(
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
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance with zero-size arrays for non-WKB or no turbulence parameterization configurations.

```julia
TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Union{Val{:TKEScheme}, Val{:TTEScheme}},
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance with zero-initialized arrays if `state.namelists.turbulence.gw_coupling == true`, otherwise the arrays are zero-size.

# Fields

  - `dtkedt0::A`: Turbulent kinetic energy impact of unresolved gravity waves.

  - `dtapedt0::A`: Turbulent available potential energy impact of unresolved gravity waves.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `turbulence_scheme`: General turbulence parameterization configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TurbulenceWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dtkedt0::A
    shearq0::A
    shearq2::A
    buoyancyq1::A
    dtapedt0::A
    convectiveq1::A
    convectiveq2::A
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBTendencies
    (; turbulence_scheme) = namelists.turbulence
    (; wkb_mode) = namelists.wkb

    @dispatch_turbulence_scheme @dispatch_wkb_mode TurbulenceWKBTendencies(
        namelists,
        domain,
        Val(turbulence_scheme),
        Val(wkb_mode),
    )
end

function TurbulenceWKBTendencies(
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
)::TurbulenceWKBTendencies
    TurbulenceWKBTendencies((zeros(0, 0, 0) for i in 1:7)...)
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBTendencies
    (; nxx, nyy, nzz) = domain
    (; gw_coupling) = namelists.turbulence

    if gw_coupling
        TurbulenceWKBTendencies(
            (zeros(nxx, nyy, nzz) for i in 1:4)...,
            (zeros(0, 0, 0) for i in 1:3)...,
        )
    else
        TurbulenceWKBTendencies((zeros(0, 0, 0) for i in 1:7)...)
    end
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TTEScheme},
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBTendencies
    (; nxx, nyy, nzz) = domain
    (; gw_coupling) = namelists.turbulence

    if gw_coupling
        TurbulenceWKBTendencies((zeros(nxx, nyy, nzz) for i in 1:7)...)
    else
        TurbulenceWKBTendencies((zeros(0, 0, 0) for i in 1:7)...)
    end
end
