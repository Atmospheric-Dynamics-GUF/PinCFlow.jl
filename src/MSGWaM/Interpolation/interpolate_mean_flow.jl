"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::U,
)
```

Interpolate the zonal wind (``u_\\mathrm{b}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``u_\\mathrm{b}`` to the location of interest, using `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``u_\\mathrm{b}`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::U,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; u) = state.variables.predictands
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (U): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr + 1 > nxx
            error(
                "Error in interpolate_mean_flow (U): ixr + 1 = ",
                ixr + 1,
                "> nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io] + dx / 2
    xl = x[ixl + io] + dx / 2

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (U): jyl = ", jyl, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow (U): jyr = ",
                jyr,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo]
    yb = y[jyb + jo]

    # Locate the closest points in vertical direction.

    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = (ztfc[ixl, jyb, kzlbd] + ztfc[ixl + 1, jyb, kzlbd]) / 2
    zlbu = (ztfc[ixl, jyb, kzlbu] + ztfc[ixl + 1, jyb, kzlbu]) / 2

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = (ztfc[ixl, jyf, kzlfd] + ztfc[ixl + 1, jyf, kzlfd]) / 2
    zlfu = (ztfc[ixl, jyf, kzlfu] + ztfc[ixl + 1, jyf, kzlfu]) / 2

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = (ztfc[ixr, jyb, kzrbd] + ztfc[ixr + 1, jyb, kzrbd]) / 2
    zrbu = (ztfc[ixr, jyb, kzrbu] + ztfc[ixr + 1, jyb, kzrbu]) / 2

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = (ztfc[ixr, jyf, kzrfd] + ztfc[ixr + 1, jyf, kzrfd]) / 2
    zrfu = (ztfc[ixr, jyf, kzrfu] + ztfc[ixr + 1, jyf, kzrfu]) / 2

    philbd = u[ixl, jyb, kzlbd]
    philbu = u[ixl, jyb, kzlbu]

    philfd = u[ixl, jyf, kzlfd]
    philfu = u[ixl, jyf, kzlfu]

    phirbd = u[ixr, jyb, kzrbd]
    phirbu = u[ixr, jyb, kzrbu]

    phirfd = u[ixr, jyf, kzrfd]
    phirfu = u[ixr, jyf, kzrfu]

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::V,
)
```

Interpolate the meridional wind (``v_\\mathrm{b}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``v_\\mathrm{b}`` to the location of interest, using `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``v_\\mathrm{b}`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::V,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; v) = state.variables.predictands
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (V): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (V): ixr = ",
                ixr,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - 1 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (V): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (V): jyf + 1 = ",
                jyf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo] + dy / 2
    yb = y[jyb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = (ztfc[ixl, jyb, kzlbd] + ztfc[ixl, jyb + 1, kzlbd]) / 2
    zlbu = (ztfc[ixl, jyb, kzlbu] + ztfc[ixl, jyb + 1, kzlbu]) / 2

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = (ztfc[ixl, jyf, kzlfd] + ztfc[ixl, jyf + 1, kzlfd]) / 2
    zlfu = (ztfc[ixl, jyf, kzlfu] + ztfc[ixl, jyf + 1, kzlfu]) / 2

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = (ztfc[ixr, jyb, kzrbd] + ztfc[ixr, jyb + 1, kzrbd]) / 2
    zrbu = (ztfc[ixr, jyb, kzrbu] + ztfc[ixr, jyb + 1, kzrbu]) / 2

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = (ztfc[ixr, jyf, kzrfd] + ztfc[ixr, jyf + 1, kzrfd]) / 2
    zrfu = (ztfc[ixr, jyf, kzrfu] + ztfc[ixr, jyf + 1, kzrfu]) / 2

    # Assign the values.

    philbd = v[ixl, jyb, kzlbd]
    philbu = v[ixl, jyb, kzlbu]

    philfd = v[ixl, jyf, kzlfd]
    philfu = v[ixl, jyf, kzlfu]

    phirbd = v[ixr, jyb, kzrbd]
    phirbu = v[ixr, jyb, kzrbu]

    phirfd = v[ixr, jyf, kzrfd]
    phirfu = v[ixr, jyf, kzrfu]

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::W,
)
```

Interpolate the vertical wind (``w_\\mathrm{b}``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``w_\\mathrm{b}`` to the location of interest, using `compute_vertical_wind` and `interpolate`. At grid points beyond the vertical boundaries, the values used in the interpolation are replaced with zeros.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``w_\\mathrm{b}`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::W,
)
    (; namelists, domain, grid) = state
    (; predictands) = state.variables
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc, topography_surface) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (W): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (W): ixr = ",
                ixr,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (W): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow (W): jyf = ",
                jyf,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo]
    yb = y[jyb + jo]

    # Locate the closest points in vertical direction.

    kzlbu = get_next_half_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = ztildetfc[ixl, jyb, kzlbd]
    zlbu = ztildetfc[ixl, jyb, kzlbu]

    kzlfu = get_next_half_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = ztildetfc[ixl, jyf, kzlfd]
    zlfu = ztildetfc[ixl, jyf, kzlfu]

    kzrbu = get_next_half_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = ztildetfc[ixr, jyb, kzrbd]
    zrbu = ztildetfc[ixr, jyb, kzrbu]

    kzrfu = get_next_half_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = ztildetfc[ixr, jyf, kzrfd]
    zrfu = ztildetfc[ixr, jyf, kzrfu]

    # Assign the values.

    if zlbu < topography_surface[ixl, jyb]
        philbd = 0.0
        philbu = 0.0
    elseif zlbd < topography_surface[ixl, jyb]
        philbd = 0.0
        philbu = compute_vertical_wind(ixl, jyb, kzlbu, predictands, grid)
    else
        philbd = compute_vertical_wind(ixl, jyb, kzlbd, predictands, grid)
        philbu = compute_vertical_wind(ixl, jyb, kzlbu, predictands, grid)
    end

    if zlfu < topography_surface[ixl, jyf]
        philfd = 0.0
        philfu = 0.0
    elseif zlfd < topography_surface[ixl, jyf]
        philfd = 0.0
        philfu = compute_vertical_wind(ixl, jyf, kzlfu, predictands, grid)
    else
        philfd = compute_vertical_wind(ixl, jyf, kzlfd, predictands, grid)
        philfu = compute_vertical_wind(ixl, jyf, kzlfu, predictands, grid)
    end

    if zrbu < topography_surface[ixr, jyb]
        phirbd = 0.0
        phirbu = 0.0
    elseif zrbd < topography_surface[ixr, jyb]
        phirbd = 0.0
        phirbu = compute_vertical_wind(ixr, jyb, kzrbu, predictands, grid)
    else
        phirbd = compute_vertical_wind(ixr, jyb, kzrbd, predictands, grid)
        phirbu = compute_vertical_wind(ixr, jyb, kzrbu, predictands, grid)
    end

    if zrfu < topography_surface[ixr, jyf]
        phirfd = 0.0
        phirfu = 0.0
    elseif zrfd < topography_surface[ixr, jyf]
        phirfd = 0.0
        phirfu = compute_vertical_wind(ixr, jyf, kzrfu, predictands, grid)
    else
        phirfd = compute_vertical_wind(ixr, jyf, kzrfd, predictands, grid)
        phirfu = compute_vertical_wind(ixr, jyf, kzrfu, predictands, grid)
    end

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDX,
)
```

Interpolate the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial x`` to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``\\partial u_\\mathrm{b} / \\partial x`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDX,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    if sizex == 1
        phi = 0.0
        return phi
    else
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        if ixl - 1 < 1
            error(
                "Error in interpolate_mean_flow (DUDX): ixl - 1 = ",
                ixl - 1,
                " < 1",
            )
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (DUDX): ixr = ",
                ixr,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (DUDX): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow (DUDX): jyf = ",
                jyf,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo]
    yb = y[jyb + jo]

    # Locate the closest points in vertical direction.

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

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, (ixl, jyb, kzlbd, kzlbu), DUDX())

    (philfd, philfu) =
        compute_derivatives(state, (ixl, jyf, kzlfd, kzlfu), DUDX())

    (phirbd, phirbu) =
        compute_derivatives(state, (ixr, jyb, kzrbd, kzrbu), DUDX())

    (phirfd, phirfu) =
        compute_derivatives(state, (ixr, jyf, kzrfd, kzrfu), DUDX())

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDY,
)
```

Interpolate the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial y`` to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``\\partial u_\\mathrm{b} / \\partial y`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDY,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DUDY): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DUDY): ixr + 1 = ",
                ixr + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io] + dx / 2
    xl = x[ixl + io] + dx / 2

    # Locate the closest points in meridional direction.
    if sizey == 1
        phi = 0.0
        return phi
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - 1 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (DUDY): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (DUDY): jyf + 1 = ",
                jyf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo] + dy / 2
    yb = y[jyb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd =
        (
            ztfc[ixl, jyb, kzlbd] +
            ztfc[ixl + 1, jyb, kzlbd] +
            ztfc[ixl, jyb + 1, kzlbd] +
            ztfc[ixl + 1, jyb + 1, kzlbd]
        ) / 4
    zlbu =
        (
            ztfc[ixl, jyb, kzlbu] +
            ztfc[ixl + 1, jyb, kzlbu] +
            ztfc[ixl, jyb + 1, kzlbu] +
            ztfc[ixl + 1, jyb + 1, kzlbu]
        ) / 4

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd =
        (
            ztfc[ixl, jyf, kzlfd] +
            ztfc[ixl + 1, jyf, kzlfd] +
            ztfc[ixl, jyf + 1, kzlfd] +
            ztfc[ixl + 1, jyf + 1, kzlfd]
        ) / 4
    zlfu =
        (
            ztfc[ixl, jyf, kzlfu] +
            ztfc[ixl + 1, jyf, kzlfu] +
            ztfc[ixl, jyf + 1, kzlfu] +
            ztfc[ixl + 1, jyf + 1, kzlfu]
        ) / 4

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd =
        (
            ztfc[ixr, jyb, kzrbd] +
            ztfc[ixr + 1, jyb, kzrbd] +
            ztfc[ixr, jyb + 1, kzrbd] +
            ztfc[ixr + 1, jyb + 1, kzrbd]
        ) / 4
    zrbu =
        (
            ztfc[ixr, jyb, kzrbu] +
            ztfc[ixr + 1, jyb, kzrbu] +
            ztfc[ixr, jyb + 1, kzrbu] +
            ztfc[ixr + 1, jyb + 1, kzrbu]
        ) / 4

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd =
        (
            ztfc[ixr, jyf, kzrfd] +
            ztfc[ixr + 1, jyf, kzrfd] +
            ztfc[ixr, jyf + 1, kzrfd] +
            ztfc[ixr + 1, jyf + 1, kzrfd]
        ) / 4
    zrfu =
        (
            ztfc[ixr, jyf, kzrfu] +
            ztfc[ixr + 1, jyf, kzrfu] +
            ztfc[ixr, jyf + 1, kzrfu] +
            ztfc[ixr + 1, jyf + 1, kzrfu]
        ) / 4

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, (ixl, jyb, kzlbd, kzlbu), DUDY())

    (philfd, philfu) =
        compute_derivatives(state, (ixl, jyf, kzlfd, kzlfu), DUDY())

    (phirbd, phirbu) =
        compute_derivatives(state, (ixr, jyb, kzrbd, kzrbu), DUDY())

    (phirfd, phirfu) =
        compute_derivatives(state, (ixr, jyf, kzrfd, kzrfu), DUDY())

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDZ,
)
```

Interpolate the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial u_\\mathrm{b} / \\partial z`` to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``\\partial u_\\mathrm{b} / \\partial z`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDZ,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DUDZ): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DUDZ): ixr + 1 = ",
                ixr + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io] + dx / 2
    xl = x[ixl + io] + dx / 2

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (DUDZ): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow (DUDZ): jyf = ",
                jyf,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo]
    yb = y[jyb + jo]

    # Locate the closest points in vertical direction.

    kzlbu = get_next_half_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = (ztildetfc[ixl, jyb, kzlbd] + ztildetfc[ixl + 1, jyb, kzlbd]) / 2
    zlbu = (ztildetfc[ixl, jyb, kzlbu] + ztildetfc[ixl + 1, jyb, kzlbu]) / 2

    kzlfu = get_next_half_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = (ztildetfc[ixl, jyf, kzlfd] + ztildetfc[ixl + 1, jyf, kzlfd]) / 2
    zlfu = (ztildetfc[ixl, jyf, kzlfu] + ztildetfc[ixl + 1, jyf, kzlfu]) / 2

    kzrbu = get_next_half_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = (ztildetfc[ixr, jyb, kzrbd] + ztildetfc[ixr + 1, jyb, kzrbd]) / 2
    zrbu = (ztildetfc[ixr, jyb, kzrbu] + ztildetfc[ixr + 1, jyb, kzrbu]) / 2

    kzrfu = get_next_half_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = (ztildetfc[ixr, jyf, kzrfd] + ztildetfc[ixr + 1, jyf, kzrfd]) / 2
    zrfu = (ztildetfc[ixr, jyf, kzrfu] + ztildetfc[ixr + 1, jyf, kzrfu]) / 2

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, (ixl, jyb, kzlbd, kzlbu), DUDZ())

    (philfd, philfu) =
        compute_derivatives(state, (ixl, jyf, kzlfd, kzlfu), DUDZ())

    (phirbd, phirbu) =
        compute_derivatives(state, (ixr, jyb, kzrbd, kzrbu), DUDZ())

    (phirfd, phirfu) =
        compute_derivatives(state, (ixr, jyf, kzrfd, kzrfu), DUDZ())

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDX,
)
```

Interpolate the zonal derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial x``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x} + \\Delta \\widehat{x} / 2`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial v_\\mathrm{b} / \\partial x`` to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``\\partial v_\\mathrm{b} / \\partial x`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDX,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        phi = 0.0
        return phi
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DVDX): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr + 1 > nxx
            error(
                "Error in interpolate_mean_flow (DVDX): ixr + 1 = ",
                ixr + 1,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io] + dx / 2
    xl = x[ixl + io] + dx / 2

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - 1 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (DVDX): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf + 1 > nyy
            error(
                "Error in interpolate_mean_flow (DVDX): jyf + 1 = ",
                jyf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo] + dy / 2
    yb = y[jyb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd =
        (
            ztfc[ixl, jyb, kzlbd] +
            ztfc[ixl + 1, jyb, kzlbd] +
            ztfc[ixl, jyb + 1, kzlbd] +
            ztfc[ixl + 1, jyb + 1, kzlbd]
        ) / 4
    zlbu =
        (
            ztfc[ixl, jyb, kzlbu] +
            ztfc[ixl + 1, jyb, kzlbu] +
            ztfc[ixl, jyb + 1, kzlbu] +
            ztfc[ixl + 1, jyb + 1, kzlbu]
        ) / 4

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd =
        (
            ztfc[ixl, jyf, kzlfd] +
            ztfc[ixl + 1, jyf, kzlfd] +
            ztfc[ixl, jyf + 1, kzlfd] +
            ztfc[ixl + 1, jyf + 1, kzlfd]
        ) / 4
    zlfu =
        (
            ztfc[ixl, jyf, kzlfu] +
            ztfc[ixl + 1, jyf, kzlfu] +
            ztfc[ixl, jyf + 1, kzlfu] +
            ztfc[ixl + 1, jyf + 1, kzlfu]
        ) / 4

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd =
        (
            ztfc[ixr, jyb, kzrbd] +
            ztfc[ixr + 1, jyb, kzrbd] +
            ztfc[ixr, jyb + 1, kzrbd] +
            ztfc[ixr + 1, jyb + 1, kzrbd]
        ) / 4
    zrbu =
        (
            ztfc[ixr, jyb, kzrbu] +
            ztfc[ixr + 1, jyb, kzrbu] +
            ztfc[ixr, jyb + 1, kzrbu] +
            ztfc[ixr + 1, jyb + 1, kzrbu]
        ) / 4

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd =
        (
            ztfc[ixr, jyf, kzrfd] +
            ztfc[ixr + 1, jyf, kzrfd] +
            ztfc[ixr, jyf + 1, kzrfd] +
            ztfc[ixr + 1, jyf + 1, kzrfd]
        ) / 4
    zrfu =
        (
            ztfc[ixr, jyf, kzrfu] +
            ztfc[ixr + 1, jyf, kzrfu] +
            ztfc[ixr, jyf + 1, kzrfu] +
            ztfc[ixr + 1, jyf + 1, kzrfu]
        ) / 4

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, (ixl, jyb, kzlbd, kzlbu), DVDX())

    (philfd, philfu) =
        compute_derivatives(state, (ixl, jyf, kzlfd, kzlfu), DVDX())

    (phirbd, phirbu) =
        compute_derivatives(state, (ixr, jyb, kzrbd, kzrbu), DVDX())

    (phirfd, phirfu) =
        compute_derivatives(state, (ixr, jyf, kzrfd, kzrfu), DVDX())

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDY,
)
```

Interpolate the meridional derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial y``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y}`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial v_\\mathrm{b} / \\partial y`` to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``\\partial v_\\mathrm{b} / \\partial y`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDY,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DVDY): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (DVDY): ixr = ",
                ixr,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        phi = 0.0
        return phi
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if jyb - 1 < 1
            error(
                "Error in interpolate_mean_flow (DVDY): jyb - 1 = ",
                jyb - 1,
                " < 1",
            )
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow (DVDY): jyf = ",
                jyf,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo]
    yb = y[jyb + jo]

    # Locate the closest points in vertical direction.

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

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, (ixl, jyb, kzlbd, kzlbu), DVDY())

    (philfd, philfu) =
        compute_derivatives(state, (ixl, jyf, kzlfd, kzlfu), DVDY())

    (phirbd, phirbu) =
        compute_derivatives(state, (ixr, jyb, kzrbd, kzrbu), DVDY())

    (phirfd, phirfu) =
        compute_derivatives(state, (ixr, jyf, kzrfd, kzrfu), DVDY())

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

"""
```julia
interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDZ,
)
```

Interpolate the vertical derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial z``) to `(xlc, ylc, zlc)`, using a trilinear-interpolation algorithm.

This method first determines the two points in ``\\widehat{x}`` and ``\\widehat{y} + \\Delta \\widehat{y} / 2`` that are closest to `xlc` and `ylc`, respectively. For each of these four horizontal positions, it then determines the two points in ``z + J \\Delta \\widehat{z} / 2`` that are closest to `zlc`. The resulting eight grid points are used to interpolate ``\\partial v_\\mathrm{b} / \\partial z`` to the location of interest, using `compute_derivatives` and `interpolate`.

# Arguments

  - `xlc`: Zonal position of interest.
  - `ylc`: Meridional position of interest.
  - `zlc`: Vertical position of interest.
  - `state`: Model state.
  - `phitype`: Mean-flow quantity to interpolate.

# Returns

  - `::Float64`: Interpolated ``\\partial v_\\mathrm{b} / \\partial z`` at the location of interest.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.compute_derivatives`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate`](@ref)
"""
function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDZ,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DVDZ): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (DVDZ): ixr = ",
                ixr,
                " > nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - 1 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow: jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf + 1 > nyy
            error(
                "Error in interpolate_mean_flow: jyf + 1 = ",
                jyf + 1,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo] + dy / 2
    yb = y[jyb + jo] + dy / 2

    # Locate the closest points in vertical direction.

    kzlbu = get_next_half_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = (ztildetfc[ixl, jyb, kzlbd] + ztildetfc[ixl, jyb + 1, kzlbd]) / 2
    zlbu = (ztildetfc[ixl, jyb, kzlbu] + ztildetfc[ixl, jyb + 1, kzlbu]) / 2

    kzlfu = get_next_half_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = (ztildetfc[ixl, jyf, kzlfd] + ztildetfc[ixl, jyf + 1, kzlfd]) / 2
    zlfu = (ztildetfc[ixl, jyf, kzlfu] + ztildetfc[ixl, jyf + 1, kzlfu]) / 2

    kzrbu = get_next_half_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = (ztildetfc[ixr, jyb, kzrbd] + ztildetfc[ixr, jyb + 1, kzrbd]) / 2
    zrbu = (ztildetfc[ixr, jyb, kzrbu] + ztildetfc[ixr, jyb + 1, kzrbu]) / 2

    kzrfu = get_next_half_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = (ztildetfc[ixr, jyf, kzrfd] + ztildetfc[ixr, jyf + 1, kzrfd]) / 2
    zrfu = (ztildetfc[ixr, jyf, kzrfu] + ztildetfc[ixr, jyf + 1, kzrfu]) / 2

    # Assign the values.

    (philbd, philbu) =
        compute_derivatives(state, (ixl, jyb, kzlbd, kzlbu), DVDZ())

    (philfd, philfu) =
        compute_derivatives(state, (ixl, jyf, kzlfd, kzlfu), DVDZ())

    (phirbd, phirbu) =
        compute_derivatives(state, (ixr, jyb, kzrbd, kzrbu), DVDZ())

    (phirfd, phirfu) =
        compute_derivatives(state, (ixr, jyf, kzrfd, kzrfu), DVDZ())

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
