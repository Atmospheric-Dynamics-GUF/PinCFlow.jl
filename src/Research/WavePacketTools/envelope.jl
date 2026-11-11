# src/Examples/WavePacketTools/envelope.jl

function envelope(parameters::NamedTuple, x::Real, y::Real, z::Real)::Real
    (; rx, ry, rz, x0, y0, z0) = parameters

    envel = 1.0

    if rx != 0.0
        envel *= exp(-(x - x0)^2 / 2 / rx^2)
    end
    if ry != 0.0
        envel *= exp(-(y - y0)^2 / 2 / ry^2)
    end
    if rz != 0.0
        envel *= exp(-(z - z0)^2 / 2 / rz^2)
    end

    return envel
end
