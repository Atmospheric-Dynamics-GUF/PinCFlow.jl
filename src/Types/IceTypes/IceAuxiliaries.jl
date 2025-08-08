"""
```julia
IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Initial states of the ice variables (currently not needed).

```julia
IceAuxiliaries(icepredictands::IcePredictands)
```

Construct an `IceAuxiliaries` instance by copying the arrays in `icepredictands`.

# Fields

- `initialn::A`: Auxiliary array for the ice-crystal number concentration.
- `initialq::A`: Auxiliary array for the ice mixing ratio.
- `initialqv::A`: Auxiliary array for the water-vapor mixing ratio.

# Arguments

- `icepredictands`: Arrays for prognostic ice variables.
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
