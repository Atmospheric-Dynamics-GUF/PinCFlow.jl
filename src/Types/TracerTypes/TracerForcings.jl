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
    tracer_setup::NoTracer,
)::TracerForcings
```

Construct a `TracerForcings` instance for configurations without tracer transport.

```julia
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::AbstractTracer,
)::TracerForcings
```

Construct a `TracerForcings` instance for configurations with tracer transport.

```julia
TracerForcings(domain::Domain, test_case::AbstractTestCase)::TracerForcings
```

Construct a `TracerForcings` instance for configurations without WKB model.

```julia
TracerForcings(domain::Domain, test_case::AbstractWKBTestCase)::TracerForcings
```

Construct a `TracerForcings` instance for configurations with tracer transport and WKB model.

# Fields

  - `chiq0::A`: Leading-order tracer forcings.

# Arguments:

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.

  - `test_case`: Teset case on which the current simulation is based.

# See also:

  - [`PinCFlow.Types.TracerTypes.TracerWKBImpact`](@ref)
"""
struct TracerForcings{A <: TracerWKBImpact}
    chiq0::A
end

function TracerForcings(namelists::Namelists, domain::Domain)::TracerForcings
    (; tracer_setup) = namelists.tracer

    return TracerForcings(namelists, domain, tracer_setup)
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::NoTracer,
)::TracerForcings
    return TracerForcings(
        [TracerWKBImpact(0, 0, 0) for field in fieldnames(TracerForcings)]...,
    )
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracer_setup::AbstractTracer,
)::TracerForcings
    (; test_case) = namelists.setting

    return TracerForcings(domain, test_case)
end

function TracerForcings(
    domain::Domain,
    test_case::AbstractTestCase,
)::TracerForcings
    return TracerForcings(
        [TracerWKBImpact(0, 0, 0) for field in fieldnames(TracerForcings)]...,
    )
end

function TracerForcings(
    domain::Domain,
    test_case::AbstractWKBTestCase,
)::TracerForcings
    (; nxx, nyy, nzz) = domain

    return TracerForcings(
        [
            TracerWKBImpact(nxx, nyy, nzz) for
            field in fieldnames(TracerForcings)
        ]...,
    )
end
