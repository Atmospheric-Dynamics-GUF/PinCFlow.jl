"""
```julia
TracerForcings{A <: TracerWKBImpact}
```

Container for `TracerWKBImpact` instance with all necessary terms for the right-hand side of the tracer equation.

```julia
TracerForcings(namelists::Namelists, domain::Domain)::TracerForcings
```

Construct a `TracerForcings` instance set according to the model configuration.

```julia
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerForcings
```

Construct a `TracerForcings` instance for configurations without tracer transport.

```julia
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerForcings
```

Construct a `TracerForcings` instance for configurations with tracer transport.

```julia
TracerForcings(domain::Domain, wkb_mode::Val{:NoWKB})::TracerForcings
```

Construct a `TracerForcings` instance for configurations without WKB model.

```julia
TracerForcings(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerForcings
```

Construct a `TracerForcings` instance for configurations with tracer transport and WKB model.

# Fields

  - `chiq0::A`: Leading-order tracer forcings.

# Arguments:

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also:

  - [`PinCFlow.Types.TracerTypes.TracerWKBImpact`](@ref)
"""
struct TracerForcings{A <: TracerWKBImpact}
    chiq0::A
end

function TracerForcings(namelists::Namelists, domain::Domain)::TracerForcings
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerForcings(
        namelists,
        domain,
        Val(tracer_setup),
    )
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:NoTracer},
)::TracerForcings
    return TracerForcings(
        [TracerWKBImpact(0, 0, 0) for field in fieldnames(TracerForcings)]...,
    )
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::Val{:TracerOn},
)::TracerForcings
    (; wkb_mode) = namelists.wkb

    @dispatch_wkb_mode return TracerForcings(domain, Val(wkb_mode))
end

function TracerForcings(domain::Domain, wkb_mode::Val{:NoWKB})::TracerForcings
    return TracerForcings(
        [TracerWKBImpact(0, 0, 0) for field in fieldnames(TracerForcings)]...,
    )
end

function TracerForcings(
    domain::Domain,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::TracerForcings
    (; nxx, nyy, nzz) = domain

    return TracerForcings(
        [
            TracerWKBImpact(nxx, nyy, nzz) for
            field in fieldnames(TracerForcings)
        ]...,
    )
end
