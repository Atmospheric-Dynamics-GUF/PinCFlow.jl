"""
```julia
TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Initial states of the tracers.

```julia
TracerAuxiliaries(tracerpredictands::TracerPredictands)::TracerAuxiliaries
```

Construct a `TracerAuxiliaries` instance by copying the arrays in `tracerpredictands`.

# Fields

  - `initialtracer::A`: Initial state of a non-dimensional tracer.

# Arguments

  - `tracerpredictands`: Tracers.
"""
struct TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialtracer::A
end

function TracerAuxiliaries(
    tracerpredictands::TracerPredictands,
)::TracerAuxiliaries
    initialtracer = copy(getfield(tracerpredictands, 1))

    return TracerAuxiliaries(initialtracer)
end
