"""
```julia
interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::Union{AbstractArray{T, 3}, Abst
```

Interpolate a scalar field to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\alpha_\\mathrm{R}`` to the location of interest, using `interpolate`.

# Arguments

  - `state`: Model state.

  - `xlc`: Zonal position of interest.

  - `ylc`: Meridional position of interest.

  - `zlc`: Vertical position of interest.

  - `field`: Scalar field to be interpolated.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_scalar end

function interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
)::Union{T, Complex{T}} where {T <: Real}
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid
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

    @ivy philbd = field[il, jb, klbd]
    @ivy philbu = field[il, jb, klbu]

    @ivy philfd = field[il, jf, klfd]
    @ivy philfu = field[il, jf, klfu]

    @ivy phirbd = field[ir, jb, krbd]
    @ivy phirbu = field[ir, jb, krbu]

    @ivy phirfd = field[ir, jf, krfd]
    @ivy phirfu = field[ir, jf, krfu]

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