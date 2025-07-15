"""
```julia
TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```
"""
struct TracerAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialtracer::A
end

"""
```julia
TracerAuxiliaries(tracerpredictands::TracerPredictands)
```
"""
function TracerAuxiliaries(tracerpredictands::TracerPredictands)
    initialtracer = copy(getfield(tracerpredictands, 1))

    return TracerAuxiliaries(initialtracer)
end
