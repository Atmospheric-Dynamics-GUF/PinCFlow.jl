"""
```julia
TracerWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Tracer impact of unresolved gravity waves.

```julia 
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance by dispatching to the appropriate `tracer_setup` method.

```julia
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance with zero-size arrays for configurations without tracer transport.

```julia 
TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance by dispatching to the appropriate `wkb_mode` method.

```julia 
TracerWKBTendencies(domain::Domain, wkb_mode::Val{:NoWKB})::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance with zero-size arrays for non-WKB configurations.

```julia
TracerWKBTendencies(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBTendencies
```

Construct a `TracerWKBTendencies` instance with zero-initialized arrays.

# Fields 

  - `dchidt0::A`: Leading-order tracer impact of unresolved gravity waves.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.
"""
struct TracerWKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dchidt0::A
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
)::TracerWKBTendencies
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerWKBTendencies(
        namelists,
        domain,
        Val(tracer_setup),
    )
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:1]...)
end

function TracerWKBTendencies(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerWKBTendencies
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerWKBTendencies(domain, Val(wkb_mode))
end

function TracerWKBTendencies(
    domain::Domain,
    wkb_mode::Val{:NoWKB},
)::TracerWKBTendencies
    return TracerWKBTendencies([zeros(0, 0, 0) for i in 1:1]...)
end

function TracerWKBTendencies(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerWKBTendencies
    (; nxx, nyy, nzz) = domain

    return TracerWKBTendencies([zeros(nxx, nyy, nzz) for i in 1:1]...)
end
