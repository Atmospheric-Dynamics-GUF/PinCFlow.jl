"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::U,
)::AbstractFloat
```

Interpolate the zonal wind (``u_\\mathrm{b}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x} + \\Delta \\hat{x} / 2`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``u_\\mathrm{b}`` to the location of interest, using `interpolate`.

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::V,
)::AbstractFloat
```

Interpolate the meridional wind (``v_\\mathrm{b}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y} + \\Delta \\hat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the zonal wind (``u_\\mathrm{b}``).

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDX,
)::AbstractFloat
```

Interpolate the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial x`` to the location of interest, using `compute_derivatives` and `interpolate`.

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDY,
)::AbstractFloat
```

Interpolate the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x} + \\Delta \\hat{x} / 2`` and ``\\hat{y} + \\Delta \\hat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial y`` to the location of interest, using `compute_derivatives` and `interpolate`.

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDZ,
)::AbstractFloat
```

Interpolate the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x} + \\Delta \\hat{x} / 2`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\hat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial z`` to the location of interest, using `compute_derivatives` and `interpolate`.

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDX,
)::AbstractFloat
```

Interpolate the zonal derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial x``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x} + \\Delta \\hat{x} / 2`` and ``\\hat{y} + \\Delta \\hat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``).

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDY,
)::AbstractFloat
```

Interpolate the meridional derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial y``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``).

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDZ,
)::AbstractFloat
```

Interpolate the vertical derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial z``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y} + \\Delta \\hat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``).

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDX,
)::AbstractFloat
```

Interpolate the zonal derivative of the tracer mixing ratio (``\\partial \\chi_\\mathrm{b} / \\partial x``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x} + \\Delta \\hat{x} / 2`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial \\chi_\\mathrm{b} / \\partial x`` to the location of interest, using `compute_derivatives` and `interpolate`.

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDY,
)::AbstractFloat
```

Interpolate the meridional derivative of the tracer mixing ratio (``\\partial \\chi_\\mathrm{b} / \\partial y``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y} + \\Delta \\hat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial \\chi_\\mathrm{b} / \\partial y`` to the location of interest, using `compute_derivatives` and `interpolate`.

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDZ,
)::AbstractFloat
```

Interpolate the vertical derivative of the tracer mixing ratio (``\\partial \\chi_\\mathrm{b} / \\partial z``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\hat{x}`` and ``\\hat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\hat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial \\chi_\\mathrm{b} / \\partial z`` to the location of interest, using `compute_derivatives` and `interpolate`.


# Arguments

  - `xlc`: Zonal position of interest.

  - `ylc`: Meridional position of interest.

  - `zlc`: Vertical position of interest.

  - `state`: Model state.

  - `phitype`: Mean-flow quantity to interpolate.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
"""
function interpolate_mean_flow end

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::U,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; u) = state.variables.predictands
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Zonal index is too large: ir + 1 = ",
                ir + 1,
                "> nxx = ",
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
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error("Meridional index is too large: jf = ", jf, " > nyy = ", nyy)
        end
    end
    yf = y[jf]
    yb = y[jb]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    zlbd = (zc[il, jb, klbd] + zc[il + 1, jb, klbd]) / 2
    zlbu = (zc[il, jb, klbu] + zc[il + 1, jb, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    zlfd = (zc[il, jf, klfd] + zc[il + 1, jf, klfd]) / 2
    zlfu = (zc[il, jf, klfu] + zc[il + 1, jf, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    zrbd = (zc[ir, jb, krbd] + zc[ir + 1, jb, krbd]) / 2
    zrbu = (zc[ir, jb, krbu] + zc[ir + 1, jb, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    zrfd = (zc[ir, jf, krfd] + zc[ir + 1, jf, krfd]) / 2
    zrfu = (zc[ir, jf, krfu] + zc[ir + 1, jf, krfu]) / 2

    philbd = u[il, jb, klbd]
    philbu = u[il, jb, klbu]

    philfd = u[il, jf, klfd]
    philfu = u[il, jf, klfu]

    phirbd = u[ir, jb, krbd]
    phirbu = u[ir, jb, krbu]

    phirfd = u[ir, jf, krfd]
    phirfu = u[ir, jf, krfu]

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::V,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; v) = state.variables.predictands
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error("Zonal index is too large: ir = ", ir, " > nxx = ", nxx)
        end
    end
    xr = x[ir]
    xl = x[il]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Meridional index is too large: jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jf] + dy / 2
    yb = y[jb] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    zlbd = (zc[il, jb, klbd] + zc[il, jb + 1, klbd]) / 2
    zlbu = (zc[il, jb, klbu] + zc[il, jb + 1, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    zlfd = (zc[il, jf, klfd] + zc[il, jf + 1, klfd]) / 2
    zlfu = (zc[il, jf, klfu] + zc[il, jf + 1, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    zrbd = (zc[ir, jb, krbd] + zc[ir, jb + 1, krbd]) / 2
    zrbu = (zc[ir, jb, krbu] + zc[ir, jb + 1, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    zrfd = (zc[ir, jf, krfd] + zc[ir, jf + 1, krfd]) / 2
    zrfu = (zc[ir, jf, krfu] + zc[ir, jf + 1, krfu]) / 2

    # Assign the values.

    philbd = v[il, jb, klbd]
    philbu = v[il, jb, klbu]

    philfd = v[il, jf, klfd]
    philfu = v[il, jf, klfu]

    phirbd = v[ir, jb, krbd]
    phirbu = v[ir, jb, krbu]

    phirfd = v[ir, jf, krfd]
    phirfu = v[ir, jf, krfu]

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDX,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    if x_size == 1
        phi = 0.0
        return phi
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il - 1 < 1
            error("Zonal index is too small: il - 1 = ", il - 1, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error("Zonal index is too large: ir = ", ir, " > nxx = ", nxx)
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
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error("Meridional index is too large: jf = ", jf, " > nyy = ", nyy)
        end
    end
    yf = y[jf]
    yb = y[jb]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    zlbd = zc[il, jb, klbd]
    zlbu = zc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    zlfd = zc[il, jf, klfd]
    zlfu = zc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    zrbd = zc[ir, jb, krbd]
    zrbu = zc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    zrfd = zc[ir, jf, krfd]
    zrfu = zc[ir, jf, krfu]

    # Assign the values.

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DUDX())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DUDX())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DUDX())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DUDX())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDY,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Zonal index is too large: ir + 1 = ",
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
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Meridional index is too large: jf + 1 = ",
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
    zlbd =
        (
            zc[il, jb, klbd] +
            zc[il + 1, jb, klbd] +
            zc[il, jb + 1, klbd] +
            zc[il + 1, jb + 1, klbd]
        ) / 4
    zlbu =
        (
            zc[il, jb, klbu] +
            zc[il + 1, jb, klbu] +
            zc[il, jb + 1, klbu] +
            zc[il + 1, jb + 1, klbu]
        ) / 4

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    zlfd =
        (
            zc[il, jf, klfd] +
            zc[il + 1, jf, klfd] +
            zc[il, jf + 1, klfd] +
            zc[il + 1, jf + 1, klfd]
        ) / 4
    zlfu =
        (
            zc[il, jf, klfu] +
            zc[il + 1, jf, klfu] +
            zc[il, jf + 1, klfu] +
            zc[il + 1, jf + 1, klfu]
        ) / 4

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    zrbd =
        (
            zc[ir, jb, krbd] +
            zc[ir + 1, jb, krbd] +
            zc[ir, jb + 1, krbd] +
            zc[ir + 1, jb + 1, krbd]
        ) / 4
    zrbu =
        (
            zc[ir, jb, krbu] +
            zc[ir + 1, jb, krbu] +
            zc[ir, jb + 1, krbu] +
            zc[ir + 1, jb + 1, krbu]
        ) / 4

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    zrfd =
        (
            zc[ir, jf, krfd] +
            zc[ir + 1, jf, krfd] +
            zc[ir, jf + 1, krfd] +
            zc[ir + 1, jf + 1, krfd]
        ) / 4
    zrfu =
        (
            zc[ir, jf, krfu] +
            zc[ir + 1, jf, krfu] +
            zc[ir, jf + 1, krfu] +
            zc[ir + 1, jf + 1, krfu]
        ) / 4

    # Assign the values.

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DUDY())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DUDY())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DUDY())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DUDY())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDZ,
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
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Zonal index is too large: ir + 1 = ",
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
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error("Meridional index is too large: jf = ", jf, " > nyy = ", nyy)
        end
    end
    yf = y[jf]
    yb = y[jb]

    # Locate the closest points in vertical direction.

    klbu = get_next_half_level(il, jb, zlc, state; dkd = 1, dku = 1)
    klbd = klbu - 1
    zlbd = (zctilde[il, jb, klbd] + zctilde[il + 1, jb, klbd]) / 2
    zlbu = (zctilde[il, jb, klbu] + zctilde[il + 1, jb, klbu]) / 2

    klfu = get_next_half_level(il, jf, zlc, state; dkd = 1, dku = 1)
    klfd = klfu - 1
    zlfd = (zctilde[il, jf, klfd] + zctilde[il + 1, jf, klfd]) / 2
    zlfu = (zctilde[il, jf, klfu] + zctilde[il + 1, jf, klfu]) / 2

    krbu = get_next_half_level(ir, jb, zlc, state; dkd = 1, dku = 1)
    krbd = krbu - 1
    zrbd = (zctilde[ir, jb, krbd] + zctilde[ir + 1, jb, krbd]) / 2
    zrbu = (zctilde[ir, jb, krbu] + zctilde[ir + 1, jb, krbu]) / 2

    krfu = get_next_half_level(ir, jf, zlc, state; dkd = 1, dku = 1)
    krfd = krfu - 1
    zrfd = (zctilde[ir, jf, krfd] + zctilde[ir + 1, jf, krfd]) / 2
    zrfu = (zctilde[ir, jf, krfu] + zctilde[ir + 1, jf, krfu]) / 2

    # Assign the values.

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DUDZ())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DUDZ())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DUDZ())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DUDZ())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDX,
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
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Zonal index is too large: ir + 1 = ",
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
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Merdional index is too large: jf + 1 = ",
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
    zlbd =
        (
            zc[il, jb, klbd] +
            zc[il + 1, jb, klbd] +
            zc[il, jb + 1, klbd] +
            zc[il + 1, jb + 1, klbd]
        ) / 4
    zlbu =
        (
            zc[il, jb, klbu] +
            zc[il + 1, jb, klbu] +
            zc[il, jb + 1, klbu] +
            zc[il + 1, jb + 1, klbu]
        ) / 4

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    zlfd =
        (
            zc[il, jf, klfd] +
            zc[il + 1, jf, klfd] +
            zc[il, jf + 1, klfd] +
            zc[il + 1, jf + 1, klfd]
        ) / 4
    zlfu =
        (
            zc[il, jf, klfu] +
            zc[il + 1, jf, klfu] +
            zc[il, jf + 1, klfu] +
            zc[il + 1, jf + 1, klfu]
        ) / 4

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    zrbd =
        (
            zc[ir, jb, krbd] +
            zc[ir + 1, jb, krbd] +
            zc[ir, jb + 1, krbd] +
            zc[ir + 1, jb + 1, krbd]
        ) / 4
    zrbu =
        (
            zc[ir, jb, krbu] +
            zc[ir + 1, jb, krbu] +
            zc[ir, jb + 1, krbu] +
            zc[ir + 1, jb + 1, krbu]
        ) / 4

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    zrfd =
        (
            zc[ir, jf, krfd] +
            zc[ir + 1, jf, krfd] +
            zc[ir, jf + 1, krfd] +
            zc[ir + 1, jf + 1, krfd]
        ) / 4
    zrfu =
        (
            zc[ir, jf, krfu] +
            zc[ir + 1, jf, krfu] +
            zc[ir, jf + 1, krfu] +
            zc[ir + 1, jf + 1, krfu]
        ) / 4

    # Assign the values.

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DVDX())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DVDX())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DVDX())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DVDX())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDY,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, zc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error("Zonal index is too large: ir = ", ir, " > nxx = ", nxx)
        end
    end
    xr = x[ir]
    xl = x[il]

    # Locate the closest points in meridional direction.
    if y_size == 1
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb - 1 < 1
            error("Meridional index is too small: jb - 1 = ", jb - 1, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error("Meridional index is too large: jf = ", jf, " > nyy = ", nyy)
        end
    end
    yf = y[jf]
    yb = y[jb]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    zlbd = zc[il, jb, klbd]
    zlbu = zc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    zlfd = zc[il, jf, klfd]
    zlfu = zc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    zrbd = zc[ir, jb, krbd]
    zrbu = zc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    zrfd = zc[ir, jf, krfd]
    zrfu = zc[ir, jf, krfu]

    # Assign the values.

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DVDY())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DVDY())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DVDY())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DVDY())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDZ,
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
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error("Zonal index is too large: ir = ", ir, " > nxx = ", nxx)
        end
    end
    xr = x[ir]
    xl = x[il]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Meridional index is too large: jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jf] + dy / 2
    yb = y[jb] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_half_level(il, jb, zlc, state; dkd = 1, dku = 1)
    klbd = klbu - 1
    zlbd = (zctilde[il, jb, klbd] + zctilde[il, jb + 1, klbd]) / 2
    zlbu = (zctilde[il, jb, klbu] + zctilde[il, jb + 1, klbu]) / 2

    klfu = get_next_half_level(il, jf, zlc, state; dkd = 1, dku = 1)
    klfd = klfu - 1
    zlfd = (zctilde[il, jf, klfd] + zctilde[il, jf + 1, klfd]) / 2
    zlfu = (zctilde[il, jf, klfu] + zctilde[il, jf + 1, klfu]) / 2

    krbu = get_next_half_level(ir, jb, zlc, state; dkd = 1, dku = 1)
    krbd = krbu - 1
    zrbd = (zctilde[ir, jb, krbd] + zctilde[ir, jb + 1, krbd]) / 2
    zrbu = (zctilde[ir, jb, krbu] + zctilde[ir, jb + 1, krbu]) / 2

    krfu = get_next_half_level(ir, jf, zlc, state; dkd = 1, dku = 1)
    krfd = krfu - 1
    zrfd = (zctilde[ir, jf, krfd] + zctilde[ir, jf + 1, krfd]) / 2
    zrfu = (zctilde[ir, jf, krfu] + zctilde[ir, jf + 1, krfu]) / 2

    # Assign the values.

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DVDZ())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DVDZ())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DVDZ())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DVDZ())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDX,
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
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Zonal index is too large: ir + 1 = ",
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
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error("Meridional index is too large: jf = ", jf, " > nyy = ", nyy)
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

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DChiDX())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DChiDX())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DChiDX())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DChiDX())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDY,
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
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error("Zonal index is too large: ir = ", ir, " > nxx = ", nxx)
        end
    end
    xr = x[ir]
    xl = x[il]

    # Locate the closest points in meridional direction.
    if y_size == 1
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Meridional index is too large: jf + 1 = ",
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

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DChiDY())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DChiDY())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DChiDY())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DChiDY())

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

@ivy function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDZ,
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
            error("Zonal index is too small: il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error("Zonal index is too large: ir = ", ir, " > nxx = ", nxx)
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
            error("Meridional index is too small: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error("Meridional index is too large: jf = ", jf, " > nyy = ", nyy)
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

    (philbd, philbu) = compute_derivatives(state, il, jb, klbd, klbu, DChiDZ())

    (philfd, philfu) = compute_derivatives(state, il, jf, klfd, klfu, DChiDZ())

    (phirbd, phirbu) = compute_derivatives(state, ir, jb, krbd, krbu, DChiDZ())

    (phirfd, phirfu) = compute_derivatives(state, ir, jf, krfd, krfu, DChiDZ())

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
