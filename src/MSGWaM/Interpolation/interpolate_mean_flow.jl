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

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``u_\\mathrm{b}`` to the location of interest, using `interpolate`.

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

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the zonal wind (``u_\\mathrm{b}``).

```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::W,
)::AbstractFloat
```

Interpolate the vertical wind (``w_\\mathrm{b}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm, and return the result.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``w_\\mathrm{b}`` to the location of interest, using `compute_vertical_wind` and `interpolate`. At grid points beyond the vertical boundaries, the values used in the interpolation are replaced with zeros.

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

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial x`` to the location of interest, using `compute_derivatives` and `interpolate`.

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

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial y`` to the location of interest, using `compute_derivatives` and `interpolate`.

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

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial z`` to the location of interest, using `compute_derivatives` and `interpolate`.

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

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``).

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

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``).

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

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. The steps that follow are analogous to those in the method for the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``).

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

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial \\chi_\\mathrm{b} / \\partial x`` to the location of interest, using `compute_derivatives` and `interpolate`.

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

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial \\chi_\\mathrm{b} / \\partial y`` to the location of interest, using `compute_derivatives` and `interpolate`.

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

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial \\chi_\\mathrm{b} / \\partial z`` to the location of interest, using `compute_derivatives` and `interpolate`.


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

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
"""
function interpolate_mean_flow end

function interpolate_mean_flow(
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
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Error in interpolate_mean_flow (U): il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Error in interpolate_mean_flow (U): ir + 1 = ",
                ir + 1,
                "> nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io] + dx / 2
    @ivy xl = x[il + io] + dx / 2

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (U): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (U): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    @ivy zlbd = (ztfc[il, jb, klbd] + ztfc[il + 1, jb, klbd]) / 2
    @ivy zlbu = (ztfc[il, jb, klbu] + ztfc[il + 1, jb, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    @ivy zlfd = (ztfc[il, jf, klfd] + ztfc[il + 1, jf, klfd]) / 2
    @ivy zlfu = (ztfc[il, jf, klfu] + ztfc[il + 1, jf, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    @ivy zrbd = (ztfc[ir, jb, krbd] + ztfc[ir + 1, jb, krbd]) / 2
    @ivy zrbu = (ztfc[ir, jb, krbu] + ztfc[ir + 1, jb, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    @ivy zrfd = (ztfc[ir, jf, krfd] + ztfc[ir + 1, jf, krfd]) / 2
    @ivy zrfu = (ztfc[ir, jf, krfu] + ztfc[ir + 1, jf, krfu]) / 2

    @ivy philbd = u[il, jb, klbd]
    @ivy philbu = u[il, jb, klbu]

    @ivy philfd = u[il, jf, klfd]
    @ivy philfu = u[il, jf, klfu]

    @ivy phirbd = u[ir, jb, krbd]
    @ivy phirbu = u[ir, jb, krbu]

    @ivy phirfd = u[ir, jf, krfd]
    @ivy phirfu = u[ir, jf, krfu]

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

function interpolate_mean_flow(
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
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_mean_flow (V): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (V): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (V): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (V): jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo] + dy / 2
    @ivy yb = y[jb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    @ivy zlbd = (ztfc[il, jb, klbd] + ztfc[il, jb + 1, klbd]) / 2
    @ivy zlbu = (ztfc[il, jb, klbu] + ztfc[il, jb + 1, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    @ivy zlfd = (ztfc[il, jf, klfd] + ztfc[il, jf + 1, klfd]) / 2
    @ivy zlfu = (ztfc[il, jf, klfu] + ztfc[il, jf + 1, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    @ivy zrbd = (ztfc[ir, jb, krbd] + ztfc[ir, jb + 1, krbd]) / 2
    @ivy zrbu = (ztfc[ir, jb, krbu] + ztfc[ir, jb + 1, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    @ivy zrfd = (ztfc[ir, jf, krfd] + ztfc[ir, jf + 1, krfd]) / 2
    @ivy zrfu = (ztfc[ir, jf, krfu] + ztfc[ir, jf + 1, krfu]) / 2

    # Assign the values.

    @ivy philbd = v[il, jb, klbd]
    @ivy philbu = v[il, jb, klbu]

    @ivy philfd = v[il, jf, klfd]
    @ivy philfu = v[il, jf, klfu]

    @ivy phirbd = v[ir, jb, krbd]
    @ivy phirbu = v[ir, jb, krbu]

    @ivy phirfd = v[ir, jf, krfd]
    @ivy phirfu = v[ir, jf, krfu]

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::W,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; predictands) = state.variables
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc, topography_surface) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_mean_flow (W): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (W): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (W): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (W): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_half_level(il, jb, zlc, state; dkd = 1)
    klbd = klbu - 1
    @ivy zlbd = ztildetfc[il, jb, klbd]
    @ivy zlbu = ztildetfc[il, jb, klbu]

    klfu = get_next_half_level(il, jf, zlc, state; dkd = 1)
    klfd = klfu - 1
    @ivy zlfd = ztildetfc[il, jf, klfd]
    @ivy zlfu = ztildetfc[il, jf, klfu]

    krbu = get_next_half_level(ir, jb, zlc, state; dkd = 1)
    krbd = krbu - 1
    @ivy zrbd = ztildetfc[ir, jb, krbd]
    @ivy zrbu = ztildetfc[ir, jb, krbu]

    krfu = get_next_half_level(ir, jf, zlc, state; dkd = 1)
    krfd = krfu - 1
    @ivy zrfd = ztildetfc[ir, jf, krfd]
    @ivy zrfu = ztildetfc[ir, jf, krfu]

    # Assign the values.

    @ivy if zlbu < topography_surface[il, jb]
        philbd = 0.0
        philbu = 0.0
    elseif zlbd < topography_surface[il, jb]
        philbd = 0.0
        philbu = compute_vertical_wind(il, jb, klbu, state)
    else
        philbd = compute_vertical_wind(il, jb, klbd, state)
        philbu = compute_vertical_wind(il, jb, klbu, state)
    end

    @ivy if zlfu < topography_surface[il, jf]
        philfd = 0.0
        philfu = 0.0
    elseif zlfd < topography_surface[il, jf]
        philfd = 0.0
        philfu = compute_vertical_wind(il, jf, klfu, state)
    else
        philfd = compute_vertical_wind(il, jf, klfd, state)
        philfu = compute_vertical_wind(il, jf, klfu, state)
    end

    @ivy if zrbu < topography_surface[ir, jb]
        phirbd = 0.0
        phirbu = 0.0
    elseif zrbd < topography_surface[ir, jb]
        phirbd = 0.0
        phirbu = compute_vertical_wind(ir, jb, krbu, state)
    else
        phirbd = compute_vertical_wind(ir, jb, krbd, state)
        phirbu = compute_vertical_wind(ir, jb, krbu, state)
    end

    @ivy if zrfu < topography_surface[ir, jf]
        phirfd = 0.0
        phirfu = 0.0
    elseif zrfd < topography_surface[ir, jf]
        phirfd = 0.0
        phirfu = compute_vertical_wind(ir, jf, krfu, state)
    else
        phirfd = compute_vertical_wind(ir, jf, krfd, state)
        phirfu = compute_vertical_wind(ir, jf, krfu, state)
    end

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDX,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    if x_size == 1
        phi = 0.0
        return phi
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il - 1 < 1
            error(
                "Error in interpolate_mean_flow (DUDX): il - 1 = ",
                il - 1,
                " < 1",
            )
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (DUDX): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DUDX): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (DUDX): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = ztfc[il, jb, klbd]
    @ivy zlbu = ztfc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = ztfc[il, jf, klfd]
    @ivy zlfu = ztfc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = ztfc[ir, jb, krbd]
    @ivy zrbu = ztfc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = ztfc[ir, jf, krfd]
    @ivy zrfu = ztfc[ir, jf, krfu]

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDY,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Error in interpolate_mean_flow (DUDY): il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DUDY): ir + 1 = ",
                ir + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io] + dx / 2
    @ivy xl = x[il + io] + dx / 2

    # Locate the closest points in meridional direction.
    if y_size == 1
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DUDY): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (DUDY): jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo] + dy / 2
    @ivy yb = y[jb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    @ivy zlbd =
        (
            ztfc[il, jb, klbd] +
            ztfc[il + 1, jb, klbd] +
            ztfc[il, jb + 1, klbd] +
            ztfc[il + 1, jb + 1, klbd]
        ) / 4
    @ivy zlbu =
        (
            ztfc[il, jb, klbu] +
            ztfc[il + 1, jb, klbu] +
            ztfc[il, jb + 1, klbu] +
            ztfc[il + 1, jb + 1, klbu]
        ) / 4

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    @ivy zlfd =
        (
            ztfc[il, jf, klfd] +
            ztfc[il + 1, jf, klfd] +
            ztfc[il, jf + 1, klfd] +
            ztfc[il + 1, jf + 1, klfd]
        ) / 4
    @ivy zlfu =
        (
            ztfc[il, jf, klfu] +
            ztfc[il + 1, jf, klfu] +
            ztfc[il, jf + 1, klfu] +
            ztfc[il + 1, jf + 1, klfu]
        ) / 4

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    @ivy zrbd =
        (
            ztfc[ir, jb, krbd] +
            ztfc[ir + 1, jb, krbd] +
            ztfc[ir, jb + 1, krbd] +
            ztfc[ir + 1, jb + 1, krbd]
        ) / 4
    @ivy zrbu =
        (
            ztfc[ir, jb, krbu] +
            ztfc[ir + 1, jb, krbu] +
            ztfc[ir, jb + 1, krbu] +
            ztfc[ir + 1, jb + 1, krbu]
        ) / 4

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    @ivy zrfd =
        (
            ztfc[ir, jf, krfd] +
            ztfc[ir + 1, jf, krfd] +
            ztfc[ir, jf + 1, krfd] +
            ztfc[ir + 1, jf + 1, krfd]
        ) / 4
    @ivy zrfu =
        (
            ztfc[ir, jf, krfu] +
            ztfc[ir + 1, jf, krfu] +
            ztfc[ir, jf + 1, krfu] +
            ztfc[ir + 1, jf + 1, krfu]
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDZ,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Error in interpolate_mean_flow (DUDZ): il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DUDZ): ir + 1 = ",
                ir + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io] + dx / 2
    @ivy xl = x[il + io] + dx / 2

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DUDZ): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (DUDZ): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_half_level(il, jb, zlc, state; dkd = 1, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = (ztildetfc[il, jb, klbd] + ztildetfc[il + 1, jb, klbd]) / 2
    @ivy zlbu = (ztildetfc[il, jb, klbu] + ztildetfc[il + 1, jb, klbu]) / 2

    klfu = get_next_half_level(il, jf, zlc, state; dkd = 1, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = (ztildetfc[il, jf, klfd] + ztildetfc[il + 1, jf, klfd]) / 2
    @ivy zlfu = (ztildetfc[il, jf, klfu] + ztildetfc[il + 1, jf, klfu]) / 2

    krbu = get_next_half_level(ir, jb, zlc, state; dkd = 1, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = (ztildetfc[ir, jb, krbd] + ztildetfc[ir + 1, jb, krbd]) / 2
    @ivy zrbu = (ztildetfc[ir, jb, krbu] + ztildetfc[ir + 1, jb, krbu]) / 2

    krfu = get_next_half_level(ir, jf, zlc, state; dkd = 1, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = (ztildetfc[ir, jf, krfd] + ztildetfc[ir + 1, jf, krfd]) / 2
    @ivy zrfu = (ztildetfc[ir, jf, krfu] + ztildetfc[ir + 1, jf, krfu]) / 2

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDX,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        phi = 0.0
        return phi
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Error in interpolate_mean_flow (DVDX): il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DVDX): ir + 1 = ",
                ir + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io] + dx / 2
    @ivy xl = x[il + io] + dx / 2

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DVDX): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (DVDX): jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo] + dy / 2
    @ivy yb = y[jb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    @ivy zlbd =
        (
            ztfc[il, jb, klbd] +
            ztfc[il + 1, jb, klbd] +
            ztfc[il, jb + 1, klbd] +
            ztfc[il + 1, jb + 1, klbd]
        ) / 4
    @ivy zlbu =
        (
            ztfc[il, jb, klbu] +
            ztfc[il + 1, jb, klbu] +
            ztfc[il, jb + 1, klbu] +
            ztfc[il + 1, jb + 1, klbu]
        ) / 4

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    @ivy zlfd =
        (
            ztfc[il, jf, klfd] +
            ztfc[il + 1, jf, klfd] +
            ztfc[il, jf + 1, klfd] +
            ztfc[il + 1, jf + 1, klfd]
        ) / 4
    @ivy zlfu =
        (
            ztfc[il, jf, klfu] +
            ztfc[il + 1, jf, klfu] +
            ztfc[il, jf + 1, klfu] +
            ztfc[il + 1, jf + 1, klfu]
        ) / 4

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    @ivy zrbd =
        (
            ztfc[ir, jb, krbd] +
            ztfc[ir + 1, jb, krbd] +
            ztfc[ir, jb + 1, krbd] +
            ztfc[ir + 1, jb + 1, krbd]
        ) / 4
    @ivy zrbu =
        (
            ztfc[ir, jb, krbu] +
            ztfc[ir + 1, jb, krbu] +
            ztfc[ir, jb + 1, krbu] +
            ztfc[ir + 1, jb + 1, krbu]
        ) / 4

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    @ivy zrfd =
        (
            ztfc[ir, jf, krfd] +
            ztfc[ir + 1, jf, krfd] +
            ztfc[ir, jf + 1, krfd] +
            ztfc[ir + 1, jf + 1, krfd]
        ) / 4
    @ivy zrfu =
        (
            ztfc[ir, jf, krfu] +
            ztfc[ir + 1, jf, krfu] +
            ztfc[ir, jf + 1, krfu] +
            ztfc[ir + 1, jf + 1, krfu]
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDY,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_mean_flow (DVDY): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (DVDY): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb - 1 < 1
            error(
                "Error in interpolate_mean_flow (DVDY): jb - 1 = ",
                jb - 1,
                " < 1",
            )
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (DVDY): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = ztfc[il, jb, klbd]
    @ivy zlbu = ztfc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = ztfc[il, jf, klfd]
    @ivy zlfu = ztfc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = ztfc[ir, jb, krbd]
    @ivy zrbu = ztfc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = ztfc[ir, jf, krfd]
    @ivy zrfu = ztfc[ir, jf, krfu]

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDZ,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_mean_flow (DVDZ): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (DVDZ): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Error in interpolate_mean_flow: jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Error in interpolate_mean_flow: jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo] + dy / 2
    @ivy yb = y[jb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_half_level(il, jb, zlc, state; dkd = 1, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = (ztildetfc[il, jb, klbd] + ztildetfc[il, jb + 1, klbd]) / 2
    @ivy zlbu = (ztildetfc[il, jb, klbu] + ztildetfc[il, jb + 1, klbu]) / 2

    klfu = get_next_half_level(il, jf, zlc, state; dkd = 1, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = (ztildetfc[il, jf, klfd] + ztildetfc[il, jf + 1, klfd]) / 2
    @ivy zlfu = (ztildetfc[il, jf, klfu] + ztildetfc[il, jf + 1, klfu]) / 2

    krbu = get_next_half_level(ir, jb, zlc, state; dkd = 1, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = (ztildetfc[ir, jb, krbd] + ztildetfc[ir, jb + 1, krbd]) / 2
    @ivy zrbu = (ztildetfc[ir, jb, krbu] + ztildetfc[ir, jb + 1, krbu]) / 2

    krfu = get_next_half_level(ir, jf, zlc, state; dkd = 1, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = (ztildetfc[ir, jf, krfd] + ztildetfc[ir, jf + 1, krfd]) / 2
    @ivy zrfu = (ztildetfc[ir, jf, krfu] + ztildetfc[ir, jf + 1, krfu]) / 2

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDX,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        phi = 0.0
        return phi
    else
        il = floor(Int, (xlc + lx / 2) / dx) + i0 - 1 - io
        if il < 1
            error("Error in interpolate_mean_flow (DChiDX): il = ", il, " < 1")
        end
        ir = il + 1
        if ir + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DChiDX): ir + 1 = ",
                ir + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io] + dx / 2
    @ivy xl = x[il + io] + dx / 2

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DChiDX): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (DChiDX): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = (ztfc[il, jb, klbd] + ztfc[il + 1, jb, klbd]) / 2
    @ivy zlbu = (ztfc[il, jb, klbu] + ztfc[il + 1, jb, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = (ztfc[il, jf, klfd] + ztfc[il + 1, jf, klfd]) / 2
    @ivy zlfu = (ztfc[il, jf, klfu] + ztfc[il + 1, jf, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = (ztfc[ir, jb, krbd] + ztfc[ir + 1, jb, krbd]) / 2
    @ivy zrbu = (ztfc[ir, jb, krbu] + ztfc[ir + 1, jb, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = (ztfc[ir, jf, krfd] + ztfc[ir + 1, jf, krfd]) / 2
    @ivy zrfu = (ztfc[ir, jf, krfu] + ztfc[ir + 1, jf, krfu]) / 2

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDY,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_mean_flow (DChiDY): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (DChiDY): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        phi = 0.0
        return phi
    else
        jb = floor(Int, (ylc + ly / 2) / dy) + j0 - 1 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DChiDY): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (DChiDY): jf + 1 = ",
                jf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo] + dy / 2
    @ivy yb = y[jb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 2, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = (ztfc[il, jf, klbd] + ztfc[il, jf + 1, klbd]) / 2
    @ivy zlbu = (ztfc[il, jf, klbu] + ztfc[il, jf + 1, klbu]) / 2

    klfu = get_next_level(il, jf, zlc, state; dkd = 2, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = (ztfc[il, jf, klfd] + ztfc[il, jf + 1, klfd]) / 2
    @ivy zlfu = (ztfc[il, jf, klfu] + ztfc[il, jf + 1, klfu]) / 2

    krbu = get_next_level(ir, jb, zlc, state; dkd = 2, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = (ztfc[ir, jb, krbd] + ztfc[ir, jb + 1, krbd]) / 2
    @ivy zrbu = (ztfc[ir, jb, krbu] + ztfc[ir, jb + 1, krbu]) / 2

    krfu = get_next_level(ir, jf, zlc, state; dkd = 2, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = (ztfc[ir, jf, krfd] + ztfc[ir, jf + 1, krfd]) / 2
    @ivy zrfu = (ztfc[ir, jf, krfu] + ztfc[ir, jf + 1, krfu]) / 2

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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDZ,
)::AbstractFloat
    (; namelists, domain, grid) = state
    (; x_size, y_size) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc) = grid

    # Locate the closest points in zonal direction.
    if x_size == 1
        il = i0
        ir = i0
    else
        il = floor(Int, (xlc + lx / 2 - dx / 2) / dx) + i0 - io
        if il < 1
            error("Error in interpolate_mean_flow (DChiDZ): il = ", il, " < 1")
        end
        ir = il + 1
        if ir > nxx
            error(
                "Error in interpolate_mean_flow (DChiDZ): ir = ",
                ir,
                " > nxx = ",
                nxx,
            )
        end
    end
    @ivy xr = x[ir + io]
    @ivy xl = x[il + io]

    # Locate the closest points in meridional direction.
    if y_size == 1
        jb = j0
        jf = j0
    else
        jb = floor(Int, (ylc + ly / 2 - dy / 2) / dy) + j0 - jo
        if jb < 1
            error("Error in interpolate_mean_flow (DChiDZ): jb = ", jb, " < 1")
        end
        jf = jb + 1
        if jf > nyy
            error(
                "Error in interpolate_mean_flow (DChiDZ): jf = ",
                jf,
                " > nyy = ",
                nyy,
            )
        end
    end
    @ivy yf = y[jf + jo]
    @ivy yb = y[jb + jo]

    # Locate the closest points in vertical direction.

    klbu = get_next_level(il, jb, zlc, state; dkd = 1, dku = 1)
    klbd = klbu - 1
    @ivy zlbd = ztildetfc[il, jb, klbd]
    @ivy zlbu = ztildetfc[il, jb, klbu]

    klfu = get_next_level(il, jf, zlc, state; dkd = 1, dku = 1)
    klfd = klfu - 1
    @ivy zlfd = ztildetfc[il, jf, klfd]
    @ivy zlfu = ztildetfc[il, jf, klfu]

    krbu = get_next_level(ir, jb, zlc, state; dkd = 1, dku = 1)
    krbd = krbu - 1
    @ivy zrbd = ztildetfc[ir, jb, krbd]
    @ivy zrbu = ztildetfc[ir, jb, krbu]

    krfu = get_next_level(ir, jf, zlc, state; dkd = 1, dku = 1)
    krfd = krfu - 1
    @ivy zrfd = ztildetfc[ir, jf, krfd]
    @ivy zrfu = ztildetfc[ir, jf, krfu]

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
