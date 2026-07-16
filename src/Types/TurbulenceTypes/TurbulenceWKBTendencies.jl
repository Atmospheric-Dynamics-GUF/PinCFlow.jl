"""
```julia
TurbulenceWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Tracer impact of unresolved gravity waves.

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
    tracer_setup::Val{:NoTracer},
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance with zero-size arrays for configurations without tracer transport.

```julia 
TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance by dispatching to the appropriate method.

```julia 
TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance with zero-size arrays for non-WKB configurations.

```julia
TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBTendencies
```

Construct a `TurbulenceWKBTendencies` instance with zero-initialized arrays if `state.namelists.tracer.leading_order_impact == true`, otherwise the arrays are zero-size.

# Fields 

  - `dchidt0::A`: Leading-order tracer impact of unresolved gravity waves.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TurbulenceWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    gw_shear::A
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceWKBTendencies
    (; turbulence_scheme) = namelists.turbulence

    @dispatch_turbulence_scheme return TurbulenceWKBTendencies(
        namelists,
        domain,
        Val(turbulence_scheme),
    )
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:NoTurbulence},
)::TurbulenceWKBTendencies
    return TurbulenceWKBTendencies([zeros(0, 0, 0) for i in 1:1]...)
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    turbulence_scheme::Val{:TKEScheme},
)::TurbulenceWKBTendencies
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TurbulenceWKBTendencies(
        namelists,
        domain,
        Val(wkb_mode),
    )
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TurbulenceWKBTendencies
    return TurbulenceWKBTendencies([zeros(0, 0, 0) for i in 1:1]...)
end

function TurbulenceWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TurbulenceWKBTendencies
    (; nxx, nyy, nzz) = domain

    return TurbulenceWKBTendencies([zeros(nxx, nyy, nzz) for i in 1:1]...)
end
