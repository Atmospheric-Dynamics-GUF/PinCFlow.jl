"""
```julia
TurbulenceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Initial states of the turbulence energies.

```julia
TurbulenceAuxiliaries(turbulencepredictands::TurbulencePredictands)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance by copying the arrays in `turbulencepredictands`.

# Fields

  - `initialturbulence::A`: Initial state of a non-dimensional turbulence energies.

# Arguments

  - `turbulencepredictands`: Turbulence energies.
"""
struct TurbulenceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialturbulence::A
end

function TurbulenceAuxiliaries(
    turbulencepredictands::TurbulencePredictands,
)::TurbulenceAuxiliaries
    return TurbulenceAuxiliaries(
        [
            copy(getfield(turbulencepredictands, field)) for
            field in 1:fieldcount(TurbulenceAuxiliaries)
        ]...,
    )
end
