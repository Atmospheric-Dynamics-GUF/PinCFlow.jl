"""
```julia
IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Initial states of the ice variables.

```julia
IceAuxiliaries(icepredictands::IcePredictands)::IceAuxiliaries
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
	iaux1::A
	iaux2::A
	iaux3::A
end

function IceAuxiliaries(icepredictands::IcePredictands)
	iaux1 = copy(getfield(icepredictands, :n))
	iaux2 = copy(getfield(icepredictands, :q))
	iaux3 = copy(getfield(icepredictands, :qv))

	return IceAuxiliaries(iaux1, iaux2, iaux3)
end
