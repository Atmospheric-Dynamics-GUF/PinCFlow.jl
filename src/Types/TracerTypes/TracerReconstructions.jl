"""
```julia
TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstruction of tracers.

The first three dimensions represent physical space, the fourth represents the physical-space dimension of the reconstruction and the fifth the two directions in which it is computed.

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
    tracer_setup::Val{:no_tracer},
)::TracerReconstructions
```

Construct a `TracerReconstructions` instance with zero-size arrays for configurations without tracer transport.

```julia
TracerReconstructions(
    domain::Domain,
    tracer_setup::Val{:tracer_on},
)::TracerReconstructions
```

Construct a `TracerReconstructions` instance with zero-initialized arrays.

# Fields

  - `chitilde::A`: Reconstructions of a non-dimensional tracer.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `tracer_setup`: General tracer-transport configuration.
"""
struct TracerReconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    chitilde::A
end

function TracerReconstructions(
    namelists::Namelists,
    domain::Domain,
)::TracerReconstructions
    (; tracer_setup) = namelists.tracer

    @dispatch_tracer_setup return TracerReconstructions(
        domain,
        Val(tracer_setup),
    )
end

function TracerReconstructions(
    domain::Domain,
    tracer_setup::Val{:no_tracer},
)::TracerReconstructions
    return TracerReconstructions(
        [
            zeros(0, 0, 0, 0, 0) for field in fieldnames(TracerReconstructions)
        ]...,
    )
end

function TracerReconstructions(
    domain::Domain,
    tracer_setup::Val{:tracer_on},
)::TracerReconstructions
    (; nxx, nyy, nzz) = domain

    return TracerReconstructions(
        [
            zeros(nxx, nyy, nzz, 3, 2) for
            field in fieldnames(TracerReconstructions)
        ]...,
    )
end
