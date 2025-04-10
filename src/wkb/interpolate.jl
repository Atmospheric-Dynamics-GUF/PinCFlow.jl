function interpolate(
    namelists::Namelists,
    philbd::AbstractFloat,
    philbu::AbstractFloat,
    philfd::AbstractFloat,
    philfu::AbstractFloat,
    phirbd::AbstractFloat,
    phirbu::AbstractFloat,
    phirfd::AbstractFloat,
    phirfu::AbstractFloat,
    zrbd::AbstractFloat,
    zrbu::AbstractFloat,
    zrfd::AbstractFloat,
    zrfu::AbstractFloat,
    zlbd::AbstractFloat,
    zlbu::AbstractFloat,
    zlfd::AbstractFloat,
    zlfu::AbstractFloat,
    zlc::AbstractFloat,
    xl::AbstractFloat,
    xr::AbstractFloat,
    xlc::AbstractFloat,
    yf::AbstractFloat,
    yb::AbstractFloat,
    ylc::AbstractFloat,
)
    (; sizex, sizey) = namelists

    # Interpolate in x.
    if sizex == 1
        phibd = philbd
        phibu = philbu

        phifd = philfd
        phifu = philfu

        zbd = zlbd
        zbu = zlbu

        zfd = zlfd
        zfu = zlfu
    else
        if xr < xl
            error("Error in interpolate: xr = ", xr, " < xl = ", xl)
        elseif xr == xl
            factor = 0.0
        elseif xlc > xr
            factor = 0.0
        elseif xlc > xl
            factor = (xr - xlc) / dx
        else
            factor = 1.0
        end

        phibd = factor * philbd + (1.0 - factor) * phirbd
        phibu = factor * philbu + (1.0 - factor) * phirbu

        phifd = factor * philfd + (1.0 - factor) * phirfd
        phifu = factor * philfu + (1.0 - factor) * phirfu

        zbd = factor * zlbd + (1.0 - factor) * zrbd
        zbu = factor * zlbu + (1.0 - factor) * zrbu

        zfd = factor * zlfd + (1.0 - factor) * zrfd
        zfu = factor * zlfu + (1.0 - factor) * zrfu
    end

    # Intepolate in y.
    if sizey == 1
        phid = phibd
        phiu = phibu

        zd = zbd
        zu = zbu
    else
        if yf < yb
            error("Error in interpolate: yf = ", yf, " < yb = ", yb)
        elseif yf == yb
            factor = 0.0
        elseif ylc > yf
            factor = 0.0
        elseif ylc > yb
            factor = (yf - ylc) / dy
        else
            factor = 1.0
        end

        phid = factor * phibd + (1.0 - factor) * phifd
        phiu = factor * phibu + (1.0 - factor) * phifu

        zd = factor * zbd + (1.0 - factor) * zfd
        zu = factor * zbu + (1.0 - factor) * zfu
    end

    # Interpolate in z.
    if zu < zd
        error("Error in interpolate: zu = ", zu, " < zd = ", zd)
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    phi = factor * phid + (1.0 - factor) * phiu

    # Return the result.
    return phi
end
