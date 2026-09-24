# src/Examples/WavePacketTools/phi.jl

"""
```julia
envelope(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the phase at ``\\left(x, y, z\\right)``.

# Arguments

  - `parameters`: Parameters of the wave-packet configuration.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.
"""
function phi end

function phi(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    (; k, l, m) = parameters

    return k * x + l * y + m * z
end
