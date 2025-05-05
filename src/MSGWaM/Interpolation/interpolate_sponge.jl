function interpolate_sponge(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey, nbz) = namelists.domain
    (; sizezz, io, jo, ko, i0, j0, k1) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid
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
    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    if kzlbd + ko > sizezz - nbz
        kzlbd = k1
        kzlbu = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]
    philbd = alphaunifiedsponge[ixl, jyb, kzlbd]
    philbu = alphaunifiedsponge[ixl, jyb, kzlbu]
    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd + ko > sizezz - nbz
        kzlfd = k1
        kzlfu = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]
    philfd = alphaunifiedsponge[ixl, jyf, kzlfd]
    philfu = alphaunifiedsponge[ixl, jyf, kzlfu]
    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd + ko > sizezz - nbz
        kzrbd = k1
        kzrbu = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]
    phirbd = alphaunifiedsponge[ixr, jyb, kzrbd]
    phirbu = alphaunifiedsponge[ixr, jyb, kzrbu]
    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd + ko > sizezz - nbz
        kzrfd = k1
        kzrfu = k1
    end
    zrfd = ztfc[ixr, jyf, kzrfd]
    zrfu = ztfc[ixr, jyf, kzrfu]
    phirfd = alphaunifiedsponge[ixr, jyf, kzrfd]
    phirfu = alphaunifiedsponge[ixr, jyf, kzrfu]

    # Interpolate.
    phi = interpolate(
        namelists;
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
