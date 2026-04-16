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

Construct a `TracerWKBIntegrals` instance by dispatching to the appropriate `tracer_setup` method.

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

Construct a `TracerWKBIntegrals` instance by dispatching to the appropriate `wkb_mode` method.

```julia 
TracerWKBIntegrals(domain::Domain, wkb_mode::Val{:NoWKB})::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-size arrays for non-WKB configurations.

```julia
TracerWKBIntegrals(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
```

Construct a `TracerWKBIntegrals` instance with zero-initialized arrays.

# Fields 

  - `uchi0::A`: Leading-order zonal tracer fluxes.

  - `vchi0::A`: Leading-order meridional tracer fluxes.

  - `wchi0::A`: Leading-order vertical tracer fluxes.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TracerWKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    uchi0::A
    vchi0::A
    wchi0::A
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBIntegrals
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerWKBIntegrals(namelists, domain, Val(tracer_setup))
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBIntegrals(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBIntegrals
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBIntegrals(domain, Val(wkb_mode))
end

function TracerWKBIntegrals(domain::Domain, wkb_mode::Val{:NoWKB})::TracerWKBIntegrals
    return TracerWKBIntegrals([zeros(0, 0, 0) for i in 1:3]...)
end

function TracerWKBIntegrals(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBIntegrals
    (; nxx, nyy, nzz) = domain

    return TracerWKBIntegrals([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
