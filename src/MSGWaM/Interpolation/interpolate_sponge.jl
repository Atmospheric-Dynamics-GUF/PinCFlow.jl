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
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid
    (; alphar) = state.sponge

    # Determine closest points in horizontal direction.
    if x_size > 1
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        ir = il + 1
    else
        il = i0
        ir = i0
    end
    @ivy xl = x[io + il]
    @ivy xr = x[io + ir]

    # Determine closest points in meridional direction.
    if y_size > 1
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        jf = jb + 1
    else
        jb = j0
        jf = j0
    end
    @ivy yb = y[jo + jb]
    @ivy yf = y[jo + jf]

    # Determine closest points in vertical direction and set interpolation
    # values.

    klbu = get_next_level(il, jb, zlc, state)
    klbd = klbu - 1
    @ivy zlbd = ztfc[il, jb, klbd]
    @ivy zlbu = ztfc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state)
    klfd = klfu - 1
    @ivy zlfd = ztfc[il, jf, klfd]
    @ivy zlfu = ztfc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state)
    krbd = krbu - 1
    @ivy zrbd = ztfc[ir, jb, krbd]
    @ivy zrbu = ztfc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state)
    krfd = krfu - 1
    @ivy zrfd = ztfc[ir, jf, krfd]
    @ivy zrfu = ztfc[ir, jf, krfu]

    @ivy philbd = alphar[il, jb, klbd]
    @ivy philbu = alphar[il, jb, klbu]

    @ivy philfd = alphar[il, jf, klfd]
    @ivy philfu = alphar[il, jf, klfu]

    @ivy phirbd = alphar[ir, jb, krbd]
    @ivy phirbu = alphar[ir, jb, krbu]

    @ivy phirfd = alphar[ir, jf, krfd]
    @ivy phirfu = alphar[ir, jf, krfu]

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
