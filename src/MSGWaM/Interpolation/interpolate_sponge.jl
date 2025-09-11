"""
```julia
interpolate_sponge(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)::AbstractFloat
```

Interpolate the Rayleigh-damping coefficient of the LHS sponge (``\\alpha_\\mathrm{R}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

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
    (; alphar) = state.sponge

    # Dermine closest points in horizontal direction.
    if sizex > 1
        ixl = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        ixr = ixl + 1
    else
        ixl = i0
        ixr = i0
    end
    @ivy xl = x[io + ixl]
    @ivy xr = x[io + ixr]

    # Determine closest points in meridional direction.
    if sizey > 1
        jyb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        jyf = jyb + 1
    else
        jyb = j0
        jyf = j0
    end
    @ivy yb = y[jo + jyb]
    @ivy yf = y[jo + jyf]

    # Determine closest points in vertical direction and set interpolation
    # values.

    kzlbu = get_next_level(ixl, jyb, zlc, state)
    kzlbd = kzlbu - 1
    @ivy zlbd = ztfc[ixl, jyb, kzlbd]
    @ivy zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, state)
    kzlfd = kzlfu - 1
    @ivy zlfd = ztfc[ixl, jyf, kzlfd]
    @ivy zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, state)
    kzrbd = kzrbu - 1
    @ivy zrbd = ztfc[ixr, jyb, kzrbd]
    @ivy zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, state)
    kzrfd = kzrfu - 1
    @ivy zrfd = ztfc[ixr, jyf, kzrfd]
    @ivy zrfu = ztfc[ixr, jyf, kzrfu]

    @ivy philbd = alphar[ixl, jyb, kzlbd]
    @ivy philbu = alphar[ixl, jyb, kzlbu]

    @ivy philfd = alphar[ixl, jyf, kzlfd]
    @ivy philfu = alphar[ixl, jyf, kzlfu]

    @ivy phirbd = alphar[ixr, jyb, kzrbd]
    @ivy phirbu = alphar[ixr, jyb, kzrbu]

    @ivy phirfd = alphar[ixr, jyf, kzrfd]
    @ivy phirfu = alphar[ixr, jyf, kzrfu]

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
