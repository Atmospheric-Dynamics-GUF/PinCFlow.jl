# examples/PinCFlowExamples/src/WavePacketTools/envelope.jl

function envelope(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    (; k, l, m, rx, ry, rz, x0, y0, z0) = parameters

    r =
        sqrt(
            (rx * k * (x - x0))^2 +
            (ry * l * (y - y0))^2 +
            (rz * m * (z - z0))^2,
        ) / pi
    if r <= 1
        return (1 + cos(pi * r)) / 2
    else
        return 0.0
    end
end
