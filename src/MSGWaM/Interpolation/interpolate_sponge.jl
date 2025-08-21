"""
```julia
interpolate_sponge(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)::AbstractFloat
```

Interpolate the Rayleigh-damping coefficient of the unified sponge (``\\alpha_\\mathrm{R}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\alpha_\\mathrm{R}`` to the location of interest, using `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.

  - `ylc`: Meridional position of interest.

  - `zlc`: Vertical position of interest.

  - `state`: Model state.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_sponge end

function interpolate_sponge(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid
    (; alphaunifiedsponge) = state.sponge

    # Dermine closest points in horizontal direction.
    if sizex > 1
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        ixr = ixl + 1
    else
        ixl = i0
        ixr = i0
    end
    xl = x[io + ixl]
    xr = x[io + ixr]

    # Determine closest points in meridional direction.
    if sizey > 1
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        jyf = jyb + 1
    else
        jyb = j0
        jyf = j0
    end
    yb = y[jo + jyb]
    yf = y[jo + jyf]

    # Determine closest points in vertical direction and set interpolation
    # values.

    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = ztfc[ixr, jyf, kzrfd]
    zrfu = ztfc[ixr, jyf, kzrfu]

    philbd = alphaunifiedsponge[ixl, jyb, kzlbd]
    philbu = alphaunifiedsponge[ixl, jyb, kzlbu]

    philfd = alphaunifiedsponge[ixl, jyf, kzlfd]
    philfu = alphaunifiedsponge[ixl, jyf, kzlfu]

    phirbd = alphaunifiedsponge[ixr, jyb, kzrbd]
    phirbu = alphaunifiedsponge[ixr, jyb, kzrbu]

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
