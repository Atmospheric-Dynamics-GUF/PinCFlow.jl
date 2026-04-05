# PinCFlowExamples.jl/src/WavePacketTools/pihat.jl

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
