function interpolate_tke end

function interpolate_tke(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere

    # Determine closest points in horizontal direction.
    if x_size > 1
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        ir = il + 1
    else
        il = i0
        ir = i0
    end
    @ivy xl = x[il]
    @ivy xr = x[ir]

    # Determine closest points in meridional direction.
    if y_size > 1
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        jf = jb + 1
    else
        jb = j0
        jf = j0
    end
    @ivy yb = y[jb]
    @ivy yf = y[jf]

    # Determine closest points in vertical direction and set interpolation
    # values.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    @ivy zlbd = zc[il, jb, klbd]
    @ivy zlbu = zc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    @ivy zlfd = zc[il, jf, klfd]
    @ivy zlfu = zc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    @ivy zrbd = zc[ir, jb, krbd]
    @ivy zrbu = zc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    @ivy zrfd = zc[ir, jf, krfd]
    @ivy zrfu = zc[ir, jf, krfu]

    @ivy philbd = tke[il, jb, klbd] / (rho[il, jb, klbd] + rhobar[il, jb, klbd])
    @ivy philbu = tke[il, jb, klbu] / (rho[il, jb, klbu] + rhobar[il, jb, klbu])

    @ivy philfd = tke[il, jf, klfd] / (rho[il, jf, klfd] + rhobar[il, jf, klfd])
    @ivy philfu = tke[il, jf, klfu] / (rho[il, jf, klfu] + rhobar[il, jf, klfu])

    @ivy phirbd = tke[ir, jb, krbd] / (rho[ir, jb, krbd] + rhobar[ir, jb, krbd])
    @ivy phirbu = tke[ir, jb, krbu] / (rho[ir, jb, krbu] + rhobar[ir, jb, krbu])

    @ivy phirfd = tke[ir, jf, krfd] / (rho[ir, jf, krfd] + rhobar[ir, jf, krfd])
    @ivy phirfu = tke[ir, jf, krfu] / (rho[ir, jf, krfu] + rhobar[ir, jf, krfu])

    # Interpolate.
    phi = interpolate(
        state;
        philbd,
        philbu,
        philfd,
        philfu,
        phirbd,
        phirbu,
        phirfd,
        phirfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlc,
        yb,
        yf,
        ylc,
        xl,
        xr,
        xlc,
    )

    return phi
end
