"""
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::U)

Interpolate zonal velocity component (U) to a given 3D location using trilinear interpolation.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::U`: Type specifier for zonal velocity component

# Returns

  - Interpolated U velocity value at the specified location
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
        if (ixl < 1)
            error("Error in interpolate_mean_flow (U): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (U): ixr = ",
                ixr,
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
        if (jyb < 1)
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::V)

Interpolate meridional velocity component (V) to a given 3D location using trilinear interpolation.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::V`: Type specifier for meridional velocity component

# Returns

  - Interpolated V velocity value at the specified location
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
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow (V): jyf = ",
                jyf,
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::W)

Interpolate vertical velocity component (W) to a given 3D location using trilinear interpolation.
Handles topography by setting velocity to zero below surface level.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::W`: Type specifier for vertical velocity component

# Returns

  - Interpolated W velocity value at the specified location
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::DUDX)

Interpolate zonal derivative of zonal velocity (∂u/∂x) to a given 3D location.
Returns zero for single-point domains in x-direction.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::DUDX`: Type specifier for ∂u/∂x derivative

# Returns

  - Interpolated ∂u/∂x value at the specified location
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::DUDY)

Interpolate meridional derivative of zonal velocity (∂u/∂y) to a given 3D location.
Returns zero for single-point domains in y-direction.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::DUDY`: Type specifier for ∂u/∂y derivative

# Returns

  - Interpolated ∂u/∂y value at the specified location
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
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (DUDY): ixr = ",
                ixr,
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::DUDZ)

Interpolate vertical derivative of zonal velocity (∂u/∂z) to a given 3D location.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::DUDZ`: Type specifier for ∂u/∂z derivative

# Returns

  - Interpolated ∂u/∂z value at the specified location
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
        if ixr > nxx
            error(
                "Error in interpolate_mean_flow (DUDZ): ixr = ",
                ixr,
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::DVDX)

Interpolate zonal derivative of meridional velocity (∂v/∂x) to a given 3D location.
Returns zero for single-point domains in x-direction.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::DVDX`: Type specifier for ∂v/∂x derivative

# Returns

  - Interpolated ∂v/∂x value at the specified location
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::DVDY)

Interpolate meridional derivative of meridional velocity (∂v/∂y) to a given 3D location.
Returns zero for single-point domains in y-direction.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::DVDY`: Type specifier for ∂v/∂y derivative

# Returns

  - Interpolated ∂v/∂y value at the specified location
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
    interpolate_mean_flow(xlc, ylc, zlc, state, phitype::DVDZ)

Interpolate vertical derivative of meridional velocity (∂v/∂z) to a given 3D location.

# Arguments

  - `xlc::AbstractFloat`: Target x-coordinate
  - `ylc::AbstractFloat`: Target y-coordinate
  - `zlc::AbstractFloat`: Target z-coordinate
  - `state::State`: Model state containing variables and grid information
  - `phitype::DVDZ`: Type specifier for ∂v/∂z derivative

# Returns

  - Interpolated ∂v/∂z value at the specified location
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
        if jyf > nyy
            error(
                "Error in interpolate_mean_flow: jyf = ",
                jyf,
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
