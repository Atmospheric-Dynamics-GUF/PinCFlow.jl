# src/Examples/WavePacketTools/vhat.jl

"""
```julia
vhat(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the meridional-wind wave amplitude at ``\\left(x, y, z\\right)``.

Consistent with the polarization relations, the meridional-wind wave amplitude is given by

```math
v_\\mathrm{w} = \\frac{i}{m N^2} \\frac{\\hat{\\omega}^2 - N^2}{\\hat{\\omega}^2 - f^2} \\left(l \\hat{\\omega} - k f\\right) b_\\mathrm{w},
```

where ``\\boldsymbol{k} = \\left(k, l, m\\right)^\\mathrm{T}`` and ``f`` are taken from `parameters` and `state.namelists.atmosphere`, respectively, and ``\\hat{\\omega}``, ``N^2`` and ``b_\\mathrm{w}`` are calculated with `omega`, `n2` and `bhat`, respectively.

# Arguments

  - `state`: Auxiliary model state.

  - `parameters`: Parameters of the wave-packet configuration.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.bhat`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.n2`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.omega`](@ref)
"""
function vhat end

function vhat(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Number
    (; coriolis_frequency) = state.namelists.atmosphere
    (; k, l, m) = parameters

    return n2(state, x, y, z) == 0.0 ? 0.0 :
           1im / m / n2(state, x, y, z) *
           (omega(state, parameters, x, y, z)^2 - n2(state, x, y, z)) /
           (omega(state, parameters, x, y, z)^2 - coriolis_frequency^2) *
           (
               l * omega(state, parameters, x, y, z) -
               1im * k * coriolis_frequency
           ) *
           bhat(state, parameters, x, y, z)
end
