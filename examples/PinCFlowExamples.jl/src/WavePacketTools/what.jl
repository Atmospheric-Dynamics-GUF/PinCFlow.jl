# examples/PinCFlowExamples.jl/src/WavePacketTools/what.jl

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
