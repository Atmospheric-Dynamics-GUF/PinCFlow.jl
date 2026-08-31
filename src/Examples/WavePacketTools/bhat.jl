# src/Examples/WavePacketTools/bhat.jl

"""
```julia
bhat(state::State, parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
```

Return the buoyancy wave amplitude at ``\\left(x, y, z\\right)``.

# Arguments

  - `state`: Auxiliary model state.

  - `parameters`: Parameters of the wave-packet configuration.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.envelope`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.n2`](@ref)
"""
function bhat end

function bhat(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
    (; a0, m) = parameters

    return a0 * n2(state, x, y, z) / m * envelope(parameters, x, y, z)
end
