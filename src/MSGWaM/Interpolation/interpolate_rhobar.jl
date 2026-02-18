function interpolate_rhobar end

function interpolate_rhobar(
    zlc::AbstractFloat,
    state::State,
    strtype::Rhobar,
)::AbstractFloat
    (; domain, grid) = state
    (; rhobar) = state.atmosphere
    (; i0, j0) = domain
    (; zc) = grid

    ku = get_next_level(i0, j0, zlc, state; dkd = 1)
    kd = ku - 1

    @ivy zd = zc[i0, j0, kd]
    @ivy zu = zc[i0, j0, ku]
    @ivy strd = rhobar[i0, j0, kd]
    @ivy stru = rhobar[i0, j0, ku]

    if zu < zd
        error(
            "Error in interpolate_rhobar (Rhobar): zu = ",
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
