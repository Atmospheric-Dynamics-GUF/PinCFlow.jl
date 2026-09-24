"""
```julia
set_p(
    model::Union{Val{:Boussinesq}, Val{:PseudoIncompressible}},
    float_type::DataType,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a zero-size array in non-compressible modes.

In these cases, the mass-weighted potential temperature is a background field: constant in Boussinesq mode, vertically varying in pseudo-incompressible mode.

```julia
set_p(
    model::Val{:Compressible},
    float_type::DataType,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a copy of ``\\bar{P} = \\bar{\\rho} \\bar{\\theta}`` in compressible mode.

In compressible mode, the mass-weighted potential temperature is a prognostic variable. Its initialization as ``P = \\bar{\\rho} \\bar{\\theta}`` means that the initial potential temperature fluctuations are such that ``\\rho \\theta = \\bar{\\rho} \\bar{\\theta}``.

# Arguments

  - `model`: Dynamic equations.

  - `float_type`: Data type of the array's elements.

  - `pbar`: Mass-weighted potential temperature.
"""
function set_p end

function set_p(
    model::Union{Val{:Boussinesq}, Val{:PseudoIncompressible}},
    float_type::DataType,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return zeros(float_type, 0, 0, 0)
end

function set_p(
    model::Val{:Compressible},
    float_type::DataType,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return copy(pbar)
end
