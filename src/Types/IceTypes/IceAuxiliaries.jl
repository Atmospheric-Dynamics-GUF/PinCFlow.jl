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
	iaux4::A
	iaux5::A
	iaux6::A
	nmax::A
	ncl::A
	clc::A

end

function IceAuxiliaries(icepredictands::IcePredictands)
	iaux1 = copy(getfield(icepredictands, :n_hom))
	iaux2 = copy(getfield(icepredictands, :q_hom))
	iaux3 = copy(getfield(icepredictands, :qv))
	iaux4 = copy(getfield(icepredictands, :n_in))
	iaux5 = copy(getfield(icepredictands, :n_het))
	iaux6 = copy(getfield(icepredictands, :q_het))
	nmax = zeros(size(iaux1))
	ncl = zeros(size(iaux1))
	clc = zeros(size(iaux1))

	return IceAuxiliaries(iaux1, iaux2, iaux3, iaux4, iaux5, iaux6, nmax, ncl, clc)
end
