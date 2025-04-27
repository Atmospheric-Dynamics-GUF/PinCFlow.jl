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
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
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
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
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
    (; nxx, nyy, io, jo, i0, j0, k0, k1) = domain
    (; lx, ly, dx, dy, x, y, ztildetfc) = grid

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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztildetfc[ixl, jyb, kzlbd]
    zlbu = ztildetfc[ixl, jyb, kzlbu]

    kzlfu = get_next_half_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztildetfc[ixl, jyf, kzlfd]
    zlfu = ztildetfc[ixl, jyf, kzlfu]

    kzrbu = get_next_half_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztildetfc[ixr, jyb, kzrbd]
    zrbu = ztildetfc[ixr, jyb, kzrbu]

    kzrfu = get_next_half_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
    zrfd = ztildetfc[ixr, jyf, kzrfd]
    zrfu = ztildetfc[ixr, jyf, kzrfu]

    # Assign the values.

    if zlbu < ztildetfc[ixl, jyb, k0 - 1]
        philbd = 0.0
        philbu = 0.0
    elseif zlbd < ztildetfc[ixl, jyb, k0 - 1]
        philbd = 0.0
        philbu = compute_vertical_wind(ixl, jyb, kzlbu, predictands, grid)
    else
        philbd = compute_vertical_wind(ixl, jyb, kzlbd, predictands, grid)
        philbu = compute_vertical_wind(ixl, jyb, kzlbu, predictands, grid)
    end

    if zlfu < ztildetfc[ixl, jyf, k0 - 1]
        philfd = 0.0
        philfu = 0.0
    elseif zlfd < ztildetfc[ixl, jyf, k0 - 1]
        philfd = 0.0
        philfu = compute_vertical_wind(ixl, jyf, kzlfu, predictands, grid)
    else
        philfd = compute_vertical_wind(ixl, jyf, kzlfd, predictands, grid)
        philfu = compute_vertical_wind(ixl, jyf, kzlfu, predictands, grid)
    end

    if zrbu < ztildetfc[ixr, jyb, k0 - 1]
        phirbd = 0.0
        phirbu = 0.0
    elseif zrbd < ztildetfc[ixr, jyb, k0 - 1]
        phirbd = 0.0
        phirbu = compute_vertical_wind(ixr, jyb, kzrbu, predictands, grid)
    else
        phirbd = compute_vertical_wind(ixr, jyb, kzrbd, predictands, grid)
        phirbu = compute_vertical_wind(ixr, jyb, kzrbu, predictands, grid)
    end

    if zrfu < ztildetfc[ixr, jyf, k0 - 1]
        phirfd = 0.0
        phirfu = 0.0
    elseif zrfd < ztildetfc[ixr, jyf, k0 - 1]
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDX,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDY,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DUDZ,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1 + 1
        kzlbd = k1 + 1
    end
    zlbd = ztildetfc[ixl, jyb, kzlbd]
    zlbu = ztildetfc[ixl, jyb, kzlbu]

    kzlfu = get_next_half_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1 + 1
        kzlfd = k1 + 1
    end
    zlfd = ztildetfc[ixl, jyf, kzlfd]
    zlfu = ztildetfc[ixl, jyf, kzlfu]

    kzrbu = get_next_half_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1 + 1
        kzrbd = k1 + 1
    end
    zrbd = ztildetfc[ixr, jyb, kzrbd]
    zrbu = ztildetfc[ixr, jyb, kzrbu]

    kzrfu = get_next_half_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1 + 1
        kzrfd = k1 + 1
    end
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDX,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDY,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1
        kzlbd = k1
    end
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1
        kzlfd = k1
    end
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1
        kzrbd = k1
    end
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1
        kzrfd = k1
    end
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

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DVDZ,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; nxx, nyy, io, jo, i0, j0, k0, k1) = domain
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
    if kzlbd > k1
        kzlbu = k1 + 1
        kzlbd = k1 + 1
    end
    zlbd = ztildetfc[ixl, jyb, kzlbd]
    zlbu = ztildetfc[ixl, jyb, kzlbu]

    kzlfu = get_next_half_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    if kzlfd > k1
        kzlfu = k1 + 1
        kzlfd = k1 + 1
    end
    zlfd = ztildetfc[ixl, jyf, kzlfd]
    zlfu = ztildetfc[ixl, jyf, kzlfu]

    kzrbu = get_next_half_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    if kzrbd > k1
        kzrbu = k1 + 1
        kzrbd = k1 + 1
    end
    zrbd = ztildetfc[ixr, jyb, kzrbd]
    zrbu = ztildetfc[ixr, jyb, kzrbu]

    kzrfu = get_next_half_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    if kzrfd > k1
        kzrfu = k1 + 1
        kzrfd = k1 + 1
    end
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
