# src/Examples/WavePacketTools/envelope.jl

function envelope(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    (; k, l, m, rx, ry, rz, x0, y0, z0, threedim) = parameters

    if threedim
        deltax = x - x0
        deltay = y - y0
        deltaz = z - z0

        env = exp(-deltaz^2 / 2 / rz^2)
        if rx != 0
            if abs(deltax) <= rx
                env *= cos(pi * deltax / 2 / rx)
            else
                env = 0
            end
        end
        if ry != 0
            if abs(deltay) <= ry
                env *= cos(pi * deltay / 2 / ry)
            else
                env = 0
            end
        end
        return env
    else
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
end
