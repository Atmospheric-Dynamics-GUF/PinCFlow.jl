"""
```julia
TurbulenceDiffusionCoefficients{A <: AbstractFloat}
```

Background values for turbulence variables.

```julia
TurbulenceDiffusionCoefficients(
    constants::Constants,
)::TurbulenceDiffusionCoefficients
```

Construct a `TurbulenceDiffusionCoefficients` instance with both fields set to ``t_\\mathrm{ref}^2 / \\left(10 L_\\mathrm{ref}^2\\right)``, where ``t_\\mathrm{ref}`` and ``L_\\mathrm{ref}`` are given by the properties `tref` and `lref` of `constants`, respectively.

# Fields

# Arguments

  - `constants`: Physical constants and reference values.
"""
struct TurbulenceDiffusionCoefficients{A <: AbstractArray{<:AbstractFloat, 3}}
    km::A
    kh::A
    kek::A
end

function TurbulenceDiffusionCoefficients(
    turbulencepredictands::TurbulencePredictands,
)::TurbulenceDiffusionCoefficients
    (; tke) = turbulencepredictands

    km = zeros(size(tke))
    kh = zeros(size(tke))
    kek = zeros(size(tke))

    return TurbulenceDiffusionCoefficients(km, kh, kek)
end
