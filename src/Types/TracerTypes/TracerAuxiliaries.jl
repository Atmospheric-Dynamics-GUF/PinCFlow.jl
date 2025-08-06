"""
```julia
TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

```julia
TracerAuxiliaries(tracerpredictands::TracerPredictands)
```
"""
struct TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialtracer::A
end

function TracerAuxiliaries(tracerpredictands::TracerPredictands)
    initialtracer = copy(getfield(tracerpredictands, 1))

    return TracerAuxiliaries(initialtracer)
end
