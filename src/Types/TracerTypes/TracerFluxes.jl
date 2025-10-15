"""
```julia
TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
```

Arrays for fluxes of tracers.

The first three dimensions represent physical space and the fourth dimension represents the flux direction.

```julia
TracerFluxes(namelists::Namelists, domain::Domain)::TracerFluxes
```

Construct a `TracerFluxes` instance with dimensions depending on the general tracer-transport configuration, by dispatching to the appropriate method.

```julia
TracerFluxes(domain::Domain, tracer_setup::NoTracer)::TracerFluxes
```

Construct a `TracerFluxes` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerFluxes(domain::Domain, tracer_setup::TracerOn)::TracerFluxes
```

Construct a `TracerFluxes` instance with zero-initialized arrays.

# Fields

  - `phichi::A`: Fluxes of a non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.
"""
struct TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phichi::A
end

function TracerFluxes(namelists::Namelists, domain::Domain)::TracerFluxes
    (; tracer_setup) = namelists.tracer

    return TracerFluxes(domain, tracer_setup)
end

function TracerFluxes(domain::Domain, tracer_setup::NoTracer)::TracerFluxes
    return TracerFluxes(
        [zeros(0, 0, 0, 0) for field in fieldnames(TracerFluxes)]...,
    )
end

function TracerFluxes(domain::Domain, tracer_setup::TracerOn)::TracerFluxes
    (; nxx, nyy, nzz) = domain

    return TracerFluxes(
        [zeros(nxx, nyy, nzz, 3) for field in fieldnames(TracerFluxes)]...,
    )
end
