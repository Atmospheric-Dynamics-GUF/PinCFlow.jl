"""
```julia 
TracerForcings{A <: TracerWKBImpact}
```

Container for `TracerWKBImpact` instance with all necessary terms for the right-hand side of the tracer equation. 

```julia 
TracerForcings(namelists::Namelists, domain::Domain)
``` 

Construct a `TracerForcings` instance set according to the model configuration.

```julia 
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracersetup::NoTracer,
)
```

Construct a `TracerForcings` instance for configurations without tracer transport. 

```julia 
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracersetup::AbstractTracer,
)
```

Construct a `TracerForcings` instance for configurations with tracer transport.

```julia 
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    testcase::AbstractTestCase,
)
```

Construct a `TracerForcings` instance for configurations without WKB model. 

```julia 
TracerForcings(
    namelists::Namelists,
    domain::Domain,
    testcase::AbstractWKBTestCase,
)
```

Construct a `TracerForcings` instance for configurations with tracer transport and WKB model.

# Arguments:

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters. 

  - `tracersetup`: General tracer-transport configuration.

  - `testcase`: Teset case on which the current simulation is based.

# See also:

  - [`PinCFlow.Types.TracerTypes.TracerWKBImpact`](@ref)

"""
struct TracerForcings{A <: TracerWKBImpact}
    chiq0::A
end

function TracerForcings(namelists::Namelists, domain::Domain)::TracerForcings
    (; tracersetup) = namelists.tracer

    return TracerForcings(namelists, domain, tracersetup)
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracersetup::NoTracer,
)::TracerForcings
    return TracerForcings(TracerWKBImpact(0, 0, 0))
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    tracersetup::AbstractTracer,
)::TracerForcings
    (; testcase) = namelists.setting

    return TracerForcings(namelists, domain, testcase)
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    testcase::AbstractTestCase,
)::TracerForcings
    return TracerForcings(TracerWKBImpact(0, 0, 0))
end

function TracerForcings(
    namelists::Namelists,
    domain::Domain,
    testcase::AbstractWKBTestCase,
)::TracerForcings
    (; nxx, nyy, nzz) = domain

    return TracerForcings(TracerWKBImpact(nxx, nyy, nzz))
end
