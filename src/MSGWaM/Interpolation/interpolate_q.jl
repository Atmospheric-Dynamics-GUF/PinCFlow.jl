function interpolate_q end

function interpolate_q(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    variable::Q00,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid
    (; q00) = state.wkb.integrals

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

    @ivy philbd = q00[il, jb, klbd]
    @ivy philbu = q00[il, jb, klbu]

    @ivy philfd = q00[il, jf, klfd]
    @ivy philfu = q00[il, jf, klfu]

    @ivy phirbd = q00[ir, jb, krbd]
    @ivy phirbu = q00[ir, jb, krbu]

    @ivy phirfd = q00[ir, jf, krfd]
    @ivy phirfu = q00[ir, jf, krfu]

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

function interpolate_q(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    variable::Q10,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid
    (; q10) = state.wkb.integrals

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

    @ivy philbd = q10[il, jb, klbd]
    @ivy philbu = q10[il, jb, klbu]

    @ivy philfd = q10[il, jf, klfd]
    @ivy philfu = q10[il, jf, klfu]

    @ivy phirbd = q10[ir, jb, krbd]
    @ivy phirbu = q10[ir, jb, krbu]

    @ivy phirfd = q10[ir, jf, krfd]
    @ivy phirfu = q10[ir, jf, krfu]

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

function interpolate_q(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    variable::Q20,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid
    (; q20) = state.wkb.integrals

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

    @ivy philbd = q20[il, jb, klbd]
    @ivy philbu = q20[il, jb, klbu]

    @ivy philfd = q20[il, jf, klfd]
    @ivy philfu = q20[il, jf, klfu]

    @ivy phirbd = q20[ir, jb, krbd]
    @ivy phirbu = q20[ir, jb, krbu]

    @ivy phirfd = q20[ir, jf, krfd]
    @ivy phirfu = q20[ir, jf, krfu]

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