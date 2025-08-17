"""
```julia
TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of tracers.

The first three dimensions represent physical space, the fourth dimension represents the direction in which the reconstruction was performed and the fifth dimension represents the two cell edges of the reconstruction.

```julia
TracerReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TracerReconstructions
```

Construct a `TracerReconstructions` instance with dimensions depending on the general tracer-transport configuration, by dispatching to the appropriate method.

```julia
TracerReconstructions(
    domain::Domain,
    tracersetup::NoTracer,
)::TracerReconstructions
```

Construct a `TracerReconstructions` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerReconstructions(
    domain::Domain,
    tracersetup::AbstractTracer,
)::TracerReconstructions
```

Construct a `TracerReconstructions` instance with zero-initialized arrays.

# Fields

  - `chitilde::A`: Reconstructions of a non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracersetup`: General tracer-transport configuration.
"""
struct TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    chitilde::A
end

function TracerReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TracerReconstructions
    (; tracersetup) = namelists.tracer

    return TracerReconstructions(domain, tracersetup)
end

function TracerReconstructions(
    domain::Domain,
    tracersetup::NoTracer,
)::TracerReconstructions
    chitilde = zeros(0, 0, 0, 0, 0)

    return TracerReconstructions(chitilde)
end

function TracerReconstructions(
    domain::Domain,
    tracersetup::AbstractTracer,
)::TracerReconstructions
    (; nxx, nyy, nzz) = domain

    chitilde = zeros(nxx, nyy, nzz, 3, 2)

    return TracerReconstructions(chitilde)
end
