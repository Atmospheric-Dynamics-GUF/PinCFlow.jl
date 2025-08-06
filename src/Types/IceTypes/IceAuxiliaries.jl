"""
```julia
IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

```julia
IceAuxiliaries(icepredictands::IcePredictands)
```
"""
struct IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialn::A
    initialq::A
    initialqv::A
end

function IceAuxiliaries(icepredictands::IcePredictands)
    initialn = copy(getfield(icepredictands, :n))
    initialq = copy(getfield(icepredictands, :q))
    initialqv = copy(getfield(icepredictands, :qv))

    return IceAuxiliaries(initialn, initialq, initialqv)
end
