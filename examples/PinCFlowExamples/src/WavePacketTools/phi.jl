# examples/PinCFlowExamples/src/WavePacketTools/phi.jl

function phi(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    (; k, l, m) = parameters

    return k * x + l * y + m * z
end
