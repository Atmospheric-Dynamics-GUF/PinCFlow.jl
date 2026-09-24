"""
```julia 
interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
)::AbstractFloat
```

Interpolate a scalar field given by the input array `field` to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate the field to the location of interest, using `interpolate`.

```julia 
interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
    derivative::DX,
)::AbstractFloat
```

Interpolate the zonal derivative of a scalar field given by the input array `field` to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result. 

This method first determines the two points in ``\\hat{x} + \\Delta \\hat{x} / 2`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate the zonal derivative of the field to the location of interest, using `compute_derivatives` and `interpolate`.

```julia 
interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
    derivative::DY,
)::AbstractFloat
```

Interpolate the meridional derivative of a scalar field given by the input array `field` to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result. 

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y} + \\Delta \\hat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate the meridional derivative of the field to the location of interest, using `compute_derivatives` and `interpolate`.

```julia 
interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
    derivative::DZ,
)::AbstractFloat
```

Interpolate the vertical derivative of a scalar field given by the input array `field` to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result. 

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in `z + J \\Delta \\hat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate the vertical derivative of the field to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `state`: Model state.

  - `xlc`: Zonal position of interest.

  - `ylc`: Meridional position of interest.

  - `zlc`: Vertical position of interest.

  - `derivative`: Direction of the interpolated derivative.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
"""
function interpolate_scalar end

@ivy function interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    # Determine closest points in zonal direction 
    if x_size > 1
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        ir = il + 1
    else
        il = i0
        ir = i0
    end
    xl = x[il]
    xr = x[ir]

    # Determine closest points in meridional direction 
    if y_size > 1
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        jf = jb + 1
    else
        jb = j0
        jf = j0
    end
    yb = y[jb]
    yf = y[jf]

    # Determine closest points in vertical direction and set interpolation
    # values 

    klbu = get_next_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    zlbd = zc[il, jb, klbd]
    zlbu = zc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    zlfd = zc[il, jf, klfd]
    zlfu = zc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    zrbd = zc[ir, jb, krbd]
    zrbu = zc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    zrfd = zc[ir, jf, krfd]
    zrfu = zc[ir, jf, krfu]

    philbd = field[il, jb, klbd]
    philbu = field[il, jb, klbu]

    philfd = field[il, jf, klfd]
    philfu = field[il, jf, klfu]

    phirbd = field[ir, jb, krbd]
    phirbu = field[ir, jb, krbu]

    phirfd = field[ir, jf, krfd]
    phirfu = field[ir, jf, krfu]

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

@ivy function interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
    derivative::DX,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        phi = 0.0
        return phi
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Error in interpolate_scalar (DX): il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Error in interpolate_scalar (DX): ir + 1 = ",
                ir + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ir] + dx / 2
    xl = x[il] + dx / 2

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_scalar (DX): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_scalar (DX): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jf]
    yb = y[jb]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    zlbd = (zc[il, jb, klbd] + zc[il + 1, jb, klbd]) / 2
    zlbu = (zc[il, jb, klbu] + zc[il + 1, jb, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    zlfd = (zc[il, jf, klfd] + zc[il + 1, jf, klfd]) / 2
    zlfu = (zc[il, jf, klfu] + zc[il + 1, jf, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    zrbd = (zc[ir, jb, krbd] + zc[ir + 1, jb, krbd]) / 2
    zrbu = (zc[ir, jb, krbu] + zc[ir + 1, jb, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    zrfd = (zc[ir, jf, krfd] + zc[ir + 1, jf, krfd]) / 2
    zrfu = (zc[ir, jf, krfu] + zc[ir + 1, jf, krfu]) / 2

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, il, jb, klbd, klbu, field, DX())

    (philfd, philfu) =
        compute_derivatives(state, il, jf, klfd, klfu, field, DX())

    (phirbd, phirbu) =
        compute_derivatives(state, ir, jb, krbd, krbu, field, DX())

    (phirfd, phirfu) =
        compute_derivatives(state, ir, jf, krfd, krfu, field, DX())

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

@ivy function interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
    derivative::DY,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_scalar (DY): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_scalar (DY): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir]
    xl = x[il]

    # Locate the closest points in meridional direction.
    if y_size == 1
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Error in interpolate_scalar (DY): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Error in interpolate_scalar (DY): jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jf] + dy / 2
    yb = y[jb] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    zlbd = (zc[il, jf, klbd] + zc[il, jf + 1, klbd]) / 2
    zlbu = (zc[il, jf, klbu] + zc[il, jf + 1, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    zlfd = (zc[il, jf, klfd] + zc[il, jf + 1, klfd]) / 2
    zlfu = (zc[il, jf, klfu] + zc[il, jf + 1, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    zrbd = (zc[ir, jb, krbd] + zc[ir, jb + 1, krbd]) / 2
    zrbu = (zc[ir, jb, krbu] + zc[ir, jb + 1, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    zrfd = (zc[ir, jf, krfd] + zc[ir, jf + 1, krfd]) / 2
    zrfu = (zc[ir, jf, krfu] + zc[ir, jf + 1, krfu]) / 2

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, il, jb, klbd, klbu, field, DY())

    (philfd, philfu) =
        compute_derivatives(state, il, jf, klfd, klfu, field, DY())

    (phirbd, phirbu) =
        compute_derivatives(state, ir, jb, krbd, krbu, field, DY())

    (phirfd, phirfu) =
        compute_derivatives(state, ir, jf, krfd, krfu, field, DY())

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

@ivy function interpolate_scalar(
    state::State,
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    field::AbstractArray{<:AbstractFloat, 3},
    derivative::DZ,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zctilde) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_scalar (DZ): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_scalar (DZ): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ir]
    xl = x[il]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_scalar (DZ): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_scalar (DZ): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jf]
    yb = y[jb]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1, dku = 1)
    klbd = klbu - 1
    zlbd = zctilde[il, jb, klbd]
    zlbu = zctilde[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 1, dku = 1)
    klfd = klfu - 1
    zlfd = zctilde[il, jf, klfd]
    zlfu = zctilde[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1, dku = 1)
    krbd = krbu - 1
    zrbd = zctilde[ir, jb, krbd]
    zrbu = zctilde[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1, dku = 1)
    krfd = krfu - 1
    zrfd = zctilde[ir, jf, krfd]
    zrfu = zctilde[ir, jf, krfu]

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, il, jb, klbd, klbu, field, DZ())

    (philfd, philfu) =
        compute_derivatives(state, il, jf, klfd, klfu, field, DZ())

    (phirbd, phirbu) =
        compute_derivatives(state, ir, jb, krbd, krbu, field, DZ())

    (phirfd, phirfu) =
        compute_derivatives(state, ir, jf, krfd, krfu, field, DZ())

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
