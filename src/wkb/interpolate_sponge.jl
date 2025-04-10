function interpolate_sponge(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)
    (; sizex, sizey) = state.namelists.domain
    (; domain, grid) = state
    (; io, jo, i0, j0, k1) = domain
    (; dx, dy, x, y, ztfc) = grid
    (; alphaunifiedsponge) = state.sponge

    # Dermine closest points in horizontal direction.
    if sizex > 1
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        ixr = ixl + 1
    else
        ixl = 1
        ixr = 1
    end
    xl = x[io + ixl]
    xr = x[io + ixr]

    # Determine closest points in meridional direction.
    if sizey > 1
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        jyf = jyb + 1
    else
        jyb = 1
        jyf = 1
    end
    yb = y[jo + jyb]
    yf = y[jo + jyf]

    # Determine closest points in vertical direction and set interpolation
    # values.
    kzlbu = kztfc(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    if kzlbd > k1
        kzlbd = k1
        kzlbu = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]
    alphalbd = alphaunifiedsponge[ixl, jyb, kzlbd]
    alphalbu = alphaunifiedsponge[ixl, jyb, kzlbu]
    kzlfu = kztfc(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfd = k1
        kzlfu = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]
    alphalfd = alphaunifiedsponge[ixl, jyf, kzlfd]
    alphalfu = alphaunifiedsponge[ixl, jyf, kzlfu]
    kzrbu = kztfc(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbd = k1
        kzrbu = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]
    alpharbd = alphaunifiedsponge[ixr, jyb, kzrbd]
    alpharbu = alphaunifiedsponge[ixr, jyb, kzrbu]
    kzrfu = kztfc(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfd = k1
        kzrfu = k1
    end
    zrfd = ztfc[ixr, jyf, kzrfd]
    zrfu = ztfc[ixr, jyf, kzrfu]
    alpharfd = alphaunifiedsponge[ixr, jyf, kzrfd]
    alpharfu = alphaunifiedsponge[ixr, jyf, kzrfu]

    # Interpolate.
    alpha = interpolate(
        namelists,
        alphalbd,
        alphalbu,
        alphalfd,
        alphalfu,
        alpharbd,
        alpharbu,
        alpharfd,
        alpharfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return alpha
end
