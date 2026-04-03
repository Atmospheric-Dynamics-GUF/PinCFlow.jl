# examples/PinCFlowExamples.jl/src/WavePacketTools/wave_action_density.jl

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
