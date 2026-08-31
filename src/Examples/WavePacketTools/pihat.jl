# src/Examples/WavePacketTools/pihat.jl

"""
```julia
pihat(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the Exner-pressure wave amplitude at ``\\left(x, y, z\\right)``.

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

  - [`PinCFlow.Examples.WavePacketTools.thetabar`](@ref)
"""
function pihat end

function pihat(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Number
    (; kappa, rsp) = state.constants
    (; m) = parameters

    return n2(state, x, y, z) == 0.0 ? 0.0 :
           kappa / rsp / thetabar(state, x, y, z) * 1im / m *
           (omega(state, parameters, x, y, z)^2 - n2(state, x, y, z)) /
           n2(state, x, y, z) * bhat(state, parameters, x, y, z)
end
