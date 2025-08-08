"""
```julia
TracerTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for the Runge-Kutta updates of tracers.

```julia
TracerTendencies(namelists::Namelists, domain::Domain)
```

Construct a `TracerTendencies` instance with dimensions depending on the general tracer-transport configuration, by dispatching to the appropriate method.

```julia
TracerTendencies(domain::Domain, tracersetup::NoTracer)
```

Construct a `TracerTendencies` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerTendencies(domain::Domain, tracersetup::AbstractTracer)
```

Construct a `TracerTendencies` instance with zero-initialized arrays.

# Fields

  - `dchi::A`: Runge-Kutta update of a non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `tracersetup`: General tracer-transport configuration.
"""
struct TracerTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dchi::A
end

function TracerTendencies(namelists::Namelists, domain::Domain)
    (; tracersetup) = namelists.tracer
    return TracerTendencies(domain, tracersetup)
end

function TracerTendencies(domain::Domain, tracersetup::NoTracer)
    dchi = zeros(0, 0, 0)

    return TracerTendencies(dchi)
end

function TracerTendencies(domain::Domain, tracersetup::AbstractTracer)
    (; nxx, nyy, nzz) = domain

    dchi = zeros(nxx, nyy, nzz)

    return TracerTendencies(dchi)
end
