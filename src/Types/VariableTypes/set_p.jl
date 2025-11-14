"""
```julia
set_p(
    model::Union{Boussinesq, PseudoIncompressible},
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a zero-size array in non-compressible modes.

In these cases, the mass-weighted potential temperature is a background field: constant in Boussinesq mode, vertically varying in pseudo-incompressible mode.

```julia
set_p(
    model::Compressible,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a copy of ``\\overline{P} = \\overline{\\rho} \\overline{\\theta}`` in compressible mode.

In compressible mode, the mass-weighted potential temperature is a prognostic variable. Its initialization as ``P = \\overline{\\rho} \\overline{\\theta}`` means that the initial potential temperature fluctuations are such that ``\\rho \\theta = \\overline{\\rho} \\overline{\\theta}``.

# Arguments:

  - `mode`: Dynamic equations.

  - `pbar`: Mass-weighted potential temperature.
"""
function set_p end

function set_p(
    model::Union{Boussinesq, PseudoIncompressible},
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return zeros(0, 0, 0)
end

function set_p(
    model::Compressible,
    pbar::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return copy(pbar)
end
