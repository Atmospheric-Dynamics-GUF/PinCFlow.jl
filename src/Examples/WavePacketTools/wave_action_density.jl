# src/Examples/WavePacketTools/wave_action_density.jl

"""
```julia
wave_action_density(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the wave-action density at ``\\left(x, y, z\\right)``.

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

  - [`PinCFlow.Examples.WavePacketTools.rhobar`](@ref)
"""
function wave_action_density end

function wave_action_density(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
    (; k, l, m) = parameters

    return n2(state, x, y, z) == 0.0 ? 0.0 :
           rhobar(state, x, y, z) / 2 *
           omega(state, parameters, x, y, z) *
           (k^2 + l^2 + m^2) / n2(state, x, y, z)^2 / (k^2 + l^2) *
           bhat(state, parameters, x, y, z)^2
end
