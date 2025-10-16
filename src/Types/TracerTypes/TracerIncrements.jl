"""
```julia
TracerIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of tracers.

```julia
TracerIncrements(namelists::Namelists, domain::Domain)::TracerIncrements
```

Construct a `TracerIncrements` instance with dimensions depending on the general tracer-transport configuration, by dispatching to the appropriate method.

```julia
TracerIncrements(domain::Domain, tracer_setup::NoTracer)::TracerIncrements
```

Construct a `TracerIncrements` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerIncrements(domain::Domain, tracer_setup::AbstractTracer)::TracerIncrements
```

Construct a `TracerIncrements` instance with zero-initialized arrays.

# Fields

  - `dchi::A`: Runge-Kutta update of a non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.
"""
struct TracerIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
    dchi::A
end

function TracerIncrements(
    namelists::Namelists,
    domain::Domain,
)::TracerIncrements
    (; tracer_setup) = namelists.tracer

    return TracerIncrements(domain, tracer_setup)
end

function TracerIncrements(
    domain::Domain,
    tracer_setup::NoTracer,
)::TracerIncrements
    return TracerIncrements(
        [zeros(0, 0, 0) for field in fieldnames(TracerIncrements)]...,
    )
end

function TracerIncrements(
    domain::Domain,
    tracer_setup::AbstractTracer,
)::TracerIncrements
    (; nxx, nyy, nzz) = domain

    return TracerIncrements(
        [zeros(nxx, nyy, nzz) for field in fieldnames(TracerIncrements)]...,
    )
end
