# src/Examples/WavePacketTools/envelope.jl

function envelope(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    if get(parameters, :version, 1) == 2
        (; rx, ry, rz, x0, y0, z0) = parameters

        abs(x - x0) > rx && rx != 0.0 && return 0.0
        abs(y - y0) > ry && ry != 0.0 && return 0.0

        envel = exp(-(z - z0)^2 / (2 * rz^2))
        if rx != 0.0
            envel *= cos(pi * (x - x0) / (2 * rx))
        end
        if ry != 0.0
            envel *= cos(pi * (y - y0) / (2 * ry))
        end
        return envel
    else
        (; k, l, m, rx, ry, rz, x0, y0, z0) = parameters

        r =
            sqrt(
                (rx * k * (x - x0))^2 +
                (ry * l * (y - y0))^2 +
                (rz * m * (z - z0))^2,
            ) / pi
        return r <= 1.0 ? (1 + cos(pi * r)) / 2 : 0
    end
end
