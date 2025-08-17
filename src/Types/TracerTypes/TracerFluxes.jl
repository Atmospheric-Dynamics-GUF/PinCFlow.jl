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
TracerFluxes(domain::Domain, tracersetup::NoTracer)::TracerFluxes
```

Construct a `TracerFluxes` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerFluxes(domain::Domain, tracersetup::AbstractTracer)::TracerFluxes
```

Construct a `TracerFluxes` instance with zero-initialized arrays.

# Fields

  - `phichi::A`: Fluxes of a non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracersetup`: General tracer-transport configuration.
"""
struct TracerFluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phichi::A
end

function TracerFluxes(namelists::Namelists, domain::Domain)::TracerFluxes
    (; tracersetup) = namelists.tracer

    return TracerFluxes(domain, tracersetup)
end

function TracerFluxes(domain::Domain, tracersetup::NoTracer)::TracerFluxes
    phichi = zeros(0, 0, 0, 0)

    return TracerFluxes(phichi)
end

function TracerFluxes(domain::Domain, tracersetup::AbstractTracer)::TracerFluxes
    (; nxx, nyy, nzz) = domain

    phichi = zeros(nxx, nyy, nzz, 3)

    return TracerFluxes(phichi)
end
