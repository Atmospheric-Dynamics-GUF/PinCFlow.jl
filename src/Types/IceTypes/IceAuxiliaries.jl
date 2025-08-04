"""
```julia
IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```
"""
struct IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialn::A
    initialq::A
    initialqv::A
end

"""
```julia
IceAuxiliaries(icepredictands::IcePredictands)
```
"""
function IceAuxiliaries(icepredictands::IcePredictands)
    initialn = copy(getfield(icepredictands, :n))
    initialq = copy(getfield(icepredictands, :q))
    initialqv = copy(getfield(icepredictands, :qv))

    return IceAuxiliaries(initialn, initialq, initialqv)
end
