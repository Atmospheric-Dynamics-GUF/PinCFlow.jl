# src/Examples/WavePacketTools/what.jl

"""
```julia
what(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the vertical-wind wave amplitude at ``\\left(x, y, z\\right)``.

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
function what end

function what(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Number
    return n2(state, x, y, z) == 0.0 ? 0.0 :
           1im * omega(state, parameters, x, y, z) / n2(state, x, y, z) *
           bhat(state, parameters, x, y, z)
end
