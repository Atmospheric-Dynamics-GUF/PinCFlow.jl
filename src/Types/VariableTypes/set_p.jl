"""
```julia
set_p(
    model::Union{Boussinesq, PseudoIncompressible},
    float_type::DataType,
    rhobar::AbstractArray{<:AbstractFloat, 3},
    thetabar::AbstractArray{<:AbstractFloat, 3},
    rhop::AbstractArray{<:AbstractFloat, 3},
    thetap::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return a zero-size array in non-compressible modes.

In these cases, the mass-weighted potential temperature is a background field: constant in Boussinesq mode, vertically varying in pseudo-incompressible mode.

```julia
set_p(
    model::Compressible,
    float_type::DataType,
    rhobar::AbstractArray{<:AbstractFloat, 3},
    thetabar::AbstractArray{<:AbstractFloat, 3},
    rhop::AbstractArray{<:AbstractFloat, 3},
    thetap::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
```

Return ``P = \\rho \\theta = \\left(\\overline{\\rho} + \\rho'\\right) \\left(\\overline{\\theta} + \\theta'\\right)`` in compressible mode.

In compressible mode, the mass-weighted potential temperature is a prognostic variable.

# Arguments:

  - `mode`: Dynamic equations.

  - `float_type`: Data type of the array's elements.

  - `rhobar`: Density background.

  - `thetabar`: Potential-temperature background.

  - `rhop`: Density fluctuations.

  - `thetap`: Potential-temperature fluctuations.
"""
function set_p end

function set_p(
    model::Union{Boussinesq, PseudoIncompressible},
    float_type::DataType,
    rhobar::AbstractArray{<:AbstractFloat, 3},
    thetabar::AbstractArray{<:AbstractFloat, 3},
    rhop::AbstractArray{<:AbstractFloat, 3},
    thetap::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return zeros(float_type, 0, 0, 0)
end

function set_p(
    model::Compressible,
    float_type::DataType,
    rhobar::AbstractArray{<:AbstractFloat, 3},
    thetabar::AbstractArray{<:AbstractFloat, 3},
    rhop::AbstractArray{<:AbstractFloat, 3},
    thetap::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return (rhobar .+ rhop) .* (thetabar .+ thetap)
end
