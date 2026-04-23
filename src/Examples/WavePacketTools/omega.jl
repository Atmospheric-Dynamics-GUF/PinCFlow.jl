# src/Examples/WavePacketTools/omega.jl

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
