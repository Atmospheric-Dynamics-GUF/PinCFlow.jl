"""
```julia
TurbulenceDiffusionCoefficients{A <: AbstractArray{<:AbstractFloat, 3}}
```

Composite type for eddy diffusion coefficients.

```julia
TurbulenceDiffusionCoefficients(
    turbulencepredictands::TurbulencePredictands,
)::TurbulenceDiffusionCoefficients
```

Construct a `TurbulenceDiffusionCoefficients` instance with zero-initialized arrays.

# Fields

  - `km::A`: Eddy diffusion coefficient for momentum ``K_\\mathrm{M}``.

  - `kh::A`: Eddy diffusion coefficient for heat ``K_\\mathrm{H}``.

  - `kek::A`: Turbulent diffusion coefficient ``K_{e_\\mathrm{k}}``.

# Arguments

  - `turbulencepredictands`: Turbulence prognostic variables.
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
