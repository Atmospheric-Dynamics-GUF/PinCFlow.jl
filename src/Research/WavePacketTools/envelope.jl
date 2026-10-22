# src/Examples/WavePacketTools/envelope.jl

function envelope(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    (; lx) = parameters

    return cos(2 * pi * x / lx) * 0.15 + (1 - 0.15)

    # return (3 + cos(pi * (z - lz / 2) / lz)) / 4 * 
    #     (3 + cos(pi * (x - lx) / lx)) / 4

    # r =
    #     sqrt(
    #         (rx * k * (x - x0))^2 +
    #         (ry * l * (y - y0))^2 +
    #         (rz * m * (z - z0))^2,
    #     ) / pi
    # if r <= 1
    #     return (1 + cos(pi * r)) / 2
    # else
    #     return 0.0
    # end
end
