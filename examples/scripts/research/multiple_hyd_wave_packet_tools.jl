
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

function goussian_envelope(alpha, x, y, z)
    (x0, y0, z0) = (x_c[alpha], y_c[alpha], z_c[alpha])
    (sigma_x, sigma_y, sigma_z) = (sigma_xc[alpha], sigma_yc[alpha], sigma_zc[alpha])

    r1 = (x - x0)^2 / (2 * sigma_x^2) + 
        (y - y0)^2 / (2 * sigma_y^2) + 
        (z - z0)^2 / (2 * sigma_z^2)

    if (abs(x-x0) <= 2.5 * sigma_x) && (abs(y-y0) <= 2.5 * sigma_y) && (abs(z-z0) <= 2.5 * sigma_z)
        return exp(-r1)

    else
        return 0.0
    end
    
end



function phi(alpha, x, y, z)
    return k[alpha] * x + l[alpha] * y + m[alpha] * z
end

function omega(alpha, x, y, z)
    return sqrt(
        (n2(x, y, z) * (k[alpha]^2 + l[alpha]^2) + coriolis_frequency^2 * m[alpha]^2) /
        (m[alpha]^2),
    )
end

function bhat(alpha, x, y, z)
    return a0[alpha] * n2(x, y, z) / m[alpha] * goussian_envelope(alpha, x, y, z)
end

function uhat(alpha, x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im / m[alpha] / n2(x, y, z) * (omega(alpha, x, y, z)^2 - n2(x, y, z)) /
           (omega(alpha, x, y, z)^2 - coriolis_frequency^2) *
           (k[alpha] * omega(alpha, x, y, z) + 1im * l[alpha] * coriolis_frequency) *
           bhat(alpha, x, y, z)
end

function vhat(alpha, x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im / m[alpha] / n2(x, y, z) * (omega(alpha, x, y, z)^2 - n2(x, y, z)) /
           (omega(alpha, x, y, z)^2 - coriolis_frequency^2) *
           (l[alpha] * omega(alpha, x, y, z) - 1im * k[alpha] * coriolis_frequency) *
           bhat(alpha, x, y, z)
end

function what(alpha, x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           1im * omega(alpha, x, y, z) / n2(x, y, z) * bhat(alpha, x, y, z)
end

function pihat(alpha, x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           kappa / rsp / thetabar(x, y, z) * 1im / m[alpha] *
           (omega(alpha, x, y, z)^2 - n2(x, y, z)) / n2(x, y, z) * bhat(alpha, x, y, z)
end

function wave_action_density(alpha, x, y, z)
    return n2(x, y, z) == 0.0 ? 0.0 :
           (rhobar(x, y, z) / 2 / omega(alpha, x, y, z)  /
           n2(x, y, z)^2) * bhat(alpha, x, y, z)^2
end
