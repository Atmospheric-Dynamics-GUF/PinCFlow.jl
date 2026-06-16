"""
```julia
TracerWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Integrals of gravity-wave induced tracer fluxes.

```julia 
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance by dispatching to the appropriate method.

```julia
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-size arrays for configurations without tracer transport.

```julia 
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance by dispatching to the appropriate method.

```julia 
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-size arrays for non-WKB configurations.

```julia
TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-initialized arrays if `state.namelists.tracer.leading_order_impact == true`, otherwise the arrays are zero-size.

# Fields 

  - `uchi0::A`: Leading-order zonal tracer fluxes.

  - `vchi0::A`: Leading-order meridional tracer fluxes.

  - `wchi0::A`: Leading-order vertical tracer fluxes.

  - `uchi1::B`: Next-order zonal tracer fluxes.

  - `vchi1::B`: Next-order meridional tracer fluxes.

  - `wchi1::B`: Next-order vertical tracer fluxes.

  - `qchi::C`: Turbulent tracer fluxes.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TracerWKBIntegrals{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
    uchi0::A
    vchi0::A
    wchi0::A
    uchi1::B
    vchi1::B
    wchi1::B
    qchi::C
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBIntegrals
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerWKBIntegrals(
        namelists,
        domain,
        Val(tracer_setup),
    )
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:7]...)
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBIntegrals
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBIntegrals(
        namelists,
        domain,
        Val(wkb_mode),
    )
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:7]...)
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
    (; nxx, nyy, nzz) = domain
    (; leading_order_impact, next_order_impact, turbulence_impact) =
        namelists.tracer

    if leading_order_impact
        uchi0 = zeros(nxx, nyy, nzz)
        vchi0 = zeros(nxx, nyy, nzz)
        wchi0 = zeros(nxx, nyy, nzz)
    else
        uchi0 = zeros(0, 0, 0)
        vchi0 = zeros(0, 0, 0)
        wchi0 = zeros(0, 0, 0)
    end
    if next_order_impact
        uchi1 = zeros(nxx, nyy, nzz)
        vchi1 = zeros(nxx, nyy, nzz)
        wchi1 = zeros(nxx, nyy, nzz)
    else
        uchi1 = zeros(0, 0, 0)
        vchi1 = zeros(0, 0, 0)
        wchi1 = zeros(0, 0, 0)
    end
    if turbulence_impact
        qchi = zeros(nxx, nyy, nzz)
    else
        qchi = zeros(0, 0, 0)
    end
    return TracerWKBIntegrals(uchi0, vchi0, wchi0, uchi1, vchi1, wchi1, qchi)
end
