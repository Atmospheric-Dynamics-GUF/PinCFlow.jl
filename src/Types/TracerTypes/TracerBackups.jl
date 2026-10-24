"""
```julia
TracerBackups{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for tracer backup copies.

```julia
TracerBackups(tracerpredictands::TracerPredictands)::TracerBackups
```

Initialize backup arrays sized according to the dimensions of the arrays in `tracerpredictands`.

# Fields

  - `chiold::A`: Tracer backup.

# Arguments

  - `tracerpredictands`: Predictands for which the backups are needed.
"""
struct TracerBackups{A <: AbstractArray{<:AbstractFloat, 3}}
    chiold::A
end

function TracerBackups(tracerpredictands::TracerPredictands)::TracerBackups
    TracerBackups(
        (
            zeros(size(getfield(tracerpredictands, field))) for
            field in fieldnames(TracerPredictands)
        )...,
    )
end
