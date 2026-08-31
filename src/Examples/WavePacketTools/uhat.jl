# src/Examples/WavePacketTools/uhat.jl

"""
```julia
uhat(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the zonal-wind wave amplitude at ``\\left(x, y, z\\right)``.

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
function uhat end

function uhat(
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
               k * omega(state, parameters, x, y, z) +
               1im * l * coriolis_frequency
           ) *
           bhat(state, parameters, x, y, z)
end
