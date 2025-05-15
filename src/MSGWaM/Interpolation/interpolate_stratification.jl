function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::N2,
)
    (; domain, grid) = state
    (; bvsstrattfc) = state.atmosphere
    (; i0, j0) = domain
    (; ztfc) = grid

    kzu = get_next_level(i0, j0, zlc, domain, grid)
    kzd = kzu - 1

    zd = ztfc[i0, j0, kzd]
    zu = ztfc[i0, j0, kzu]
    strd = bvsstrattfc[i0, j0, kzd]
    stru = bvsstrattfc[i0, j0, kzu]

    if zu < zd
        error(
            "Error in interpolate_stratification (N2): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    str = factor * strd + (1.0 - factor) * stru

    return str
end

function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::DN2DZ,
)
    (; domain, grid) = state
    (; bvsstrattfc) = state.atmosphere
    (; i0, j0) = domain
    (; dz, ztildetfc, jac) = grid

    kzu = get_next_half_level(i0, j0, zlc, domain, grid)
    kzd = kzu - 1

    zd = ztildetfc[i0, j0, kzd]
    zu = ztildetfc[i0, j0, kzu]

    strd =
        (bvsstrattfc[i0, j0, kzd + 1] - bvsstrattfc[i0, j0, kzd]) / (
            2.0 * jac[i0, j0, kzd] * jac[i0, j0, kzd + 1] /
            (jac[i0, j0, kzd] + jac[i0, j0, kzd + 1])
        ) / dz
    stru =
        (bvsstrattfc[i0, j0, kzu + 1] - bvsstrattfc[i0, j0, kzu]) / (
            2.0 * jac[i0, j0, kzu] * jac[i0, j0, kzu + 1] /
            (jac[i0, j0, kzu] + jac[i0, j0, kzu + 1])
        ) / dz

    if zu < zd
        error(
            "Error in interpolate_stratification (DN2DZ): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    str = factor * strd + (1.0 - factor) * stru

    return str
end
