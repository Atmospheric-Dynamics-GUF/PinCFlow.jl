"""
```julia
TurbulenceBackups{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for turbulence backup copies.

```julia
TurbulenceBackups(
    turbulencepredictands::TurbulencePredictands,
)::TurbulenceBackups
```

Initialize backup arrays sized according to the dimensions of the arrays in `turbulencepredictands`.

# Fields

  - `tkeold::A`: Backup of the turbulent kinetic energy.

# Arguments

  - `turbulencepredictands`: Predictands for which the backups are needed.
"""
struct TurbulenceBackups{A <: AbstractArray{<:AbstractFloat, 3}}
    tkeold::A
end

function TurbulenceBackups(
    turbulencepredictands::TurbulencePredictands,
)::TurbulenceBackups
    TurbulenceBackups(
        (
            zeros(size(getfield(turbulencepredictands, field))) for
            field in fieldnames(TurbulencePredictands)
        )...,
    )
end
