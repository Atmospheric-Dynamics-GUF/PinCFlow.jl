
function ijk(x, y, z)
    i = argmin(abs.(x .- auxiliary_state.grid.x .* lref))
    j = argmin(abs.(y .- auxiliary_state.grid.y .* lref))
    k = argmin(abs.(z .- auxiliary_state.grid.zc[i, j, :] .* lref))

    return CartesianIndex(i, j, k)
end

function rhobar(x, y, z)
    return auxiliary_state.atmosphere.rhobar[ijk(x, y, z)] .* rhoref
end

function thetabar(x, y, z)
    return auxiliary_state.atmosphere.thetabar[ijk(x, y, z)] .* thetaref
end

function n2(x, y, z)
    return auxiliary_state.atmosphere.n2[ijk(x, y, z)] ./ tref .^ 2
end

function envelope(x, y, z)
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

function goussian_envelope(x, y, z)
    r1 = (x - x0)^2 / (2 * sigma_x^2) + 
         (y - y0)^2 / (2 * sigma_y^2) + 
         (z - z0)^2 / (2 * sigma_z^2)

    #if (abs(x-x0) <= 2.5 * sigma_x) && (abs(y-y0) <= 2.5 * sigma_y) && (abs(z-z0) <= 2.5 * sigma_x)
    return exp(-r1)
    #else
    #    return 0.0
    #end
end


function phi(x, y, z)
    return k * x + l * y + m * z
end

function omega(x, y, z)
    return -sqrt(
        (n2(x, y, z) * (k^2 + l^2) + coriolis_frequency^2 * m^2) /
        (k^2 + l^2 + m^2),
    )
end

function bhat(x, y, z)
    return a0 * n2(x, y, z) / m * envelope(x, y, z)
end

function uhat(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im / m / n2(x, y, z) * (omega(x, y, z)^2 - n2(x, y, z)) /
           (omega(x, y, z)^2 - coriolis_frequency^2) *
           (k * omega(x, y, z) + 1im * l * coriolis_frequency) *
           bhat(x, y, z)
end

function vhat(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im / m / n2(x, y, z) * (omega(x, y, z)^2 - n2(x, y, z)) /
           (omega(x, y, z)^2 - coriolis_frequency^2) *
           (l * omega(x, y, z) - 1im * k * coriolis_frequency) *
           bhat(x, y, z)
end

function what(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im * omega(x, y, z) / n2(x, y, z) * bhat(x, y, z)
end

function pihat(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           kappa / rsp / thetabar(x, y, z) * 1im / m *
           (omega(x, y, z)^2 - n2(x, y, z)) / n2(x, y, z) * bhat(x, y, z)
end

function wave_action_density(x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           rhobar(x, y, z) / 2 * omega(x, y, z) * (k^2 + l^2 + m^2) /
           n2(x, y, z)^2 / (k^2 + l^2) * bhat(x, y, z)^2
end
