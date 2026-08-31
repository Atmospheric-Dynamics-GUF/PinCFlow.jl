# src/Examples/WavePacketTools/omega.jl

"""
```julia
omega(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
```

Return the intrinsic frequency at ``\\left(x, y, z\\right)``.

# Arguments

  - `state`: Auxiliary model state (needed by `n2`).

  - `parameters`: Parameters of the wave-packet configuration.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.n2`](@ref)
"""
function omega end

function omega(
    state::State,
    parameters::NamedTuple,
    x::Real,
    y::Real,
    z::Real,
)::Real
    (; coriolis_frequency) = state.namelists.atmosphere
    (; k, l, m) = parameters

    return -sqrt(
        (n2(state, x, y, z) * (k^2 + l^2) + coriolis_frequency^2 * m^2) /
        (k^2 + l^2 + m^2),
    )
end
