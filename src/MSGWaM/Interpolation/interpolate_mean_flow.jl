function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    flwtype::AbstractVariable,
)
    (; namelists, domain, grid) = State
    (; predictands) = state.variables
    return interpolate_mean_flow(
        xlc,
        ylc,
        zlc,
        namelists,
        domain,
        grid,
        predictands,
        flwtype,
    )
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::U,
)
    (; sizex, sizey) = namelists.domain
    (; u) = predictands
    (; io, jo, i0, j0, k1) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - io
        if (ixl < 1)
            error("Error in interpolate_mean_flow (U): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error("Error in interpolate_mean_flow (U): ixr = ", ixr, "> nxx = ", nxx)
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
            error("Error in interpolate_mean_flow (U): jyr = ", jyr, " > nyy = ", nyy)
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

    flwlbd = u[ixl, jyb, kzlbd]
    flwlbu = u[ixl, jyb, kzlbu]

    flwlfd = u[ixl, jyf, kzlfd]
    flwlfu = u[ixl, jyf, kzlfu]

    flwrbd = u[ixr, jyb, kzrbd]
    flwrbu = u[ixr, jyb, kzrbu]

    flwrfd = u[ixr, jyf, kzrfd]
    flwrfu = u[ixr, jyf, kzrfu]

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::V,
)
    (; sizex, sizey) = namelists.domain
    (; v) = predictands
    (; io, jo, i0, j0, k1) = domain
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
            error("Error in interpolate_mean_flow (V): ixr = ", ixr, " > nxx = ", nxx)
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow (V): jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error("Error in interpolate_mean_flow (V): jyf = ", jyf, " > nyy = ", nyy)
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

    flwlbd = v[ixl, jyb, kzlbd]
    flwlbu = v[ixl, jyb, kzlbu]

    flwlfd = v[ixl, jyf, kzlfd]
    flwlfu = v[ixl, jyf, kzlfu]

    flwrbd = v[ixr, jyb, kzrbd]
    flwrbu = v[ixr, jyb, kzrbu]

    flwrfd = v[ixr, jyf, kzrfd]
    flwrfu = v[ixr, jyf, kzrfu]

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::W,
)
    (; sizex, sizey) = namelists.domain
    (; io, jo, i0, j0, k0, k1) = domain
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
            error("Error in interpolate_mean_flow (W): ixr = ", ixr, " > nxx = ", nxx)
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
            error("Error in interpolate_mean_flow (W): jyf = ", jyf, " > nyy = ", nyy)
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
        flwlbd = 0.0
        flwlbu = 0.0
    elseif zlbd < ztildetfc[ixl, jyb, k0 - 1]
        flwlbd = 0.0
        flwlbu = compute_vertical_wind(ixl, jyb, kzlbu, predictands, grid)
    else
        flwlbd = compute_vertical_wind(ixl, jyb, kzlbd, predictands, grid)
        flwlbu = compute_vertical_wind(ixl, jyb, kzlbu, predictands, grid)
    end

    if zlfu < ztildetfc[ixl, jyf, k0 - 1]
        flwlfd = 0.0
        flwlfu = 0.0
    elseif zlfd < ztildetfc[ixl, jyf, k0 - 1]
        flwlfd = 0.0
        flwlfu = compute_vertical_wind(ixl, jyf, kzlfu, predictands, grid)
    else
        flwlfd = compute_vertical_wind(ixl, jyf, kzlfd, predictands, grid)
        flwlfu = compute_vertical_wind(ixl, jyf, kzlfu, predictands, grid)
    end

    if zrbu < ztildetfc[ixr, jyb, k0 - 1]
        flwrbd = 0.0
        flwrbu = 0.0
    elseif zrbd < ztildetfc[ixr, jyb, k0 - 1]
        flwrbd = 0.0
        flwrbu = compute_vertical_wind(ixr, jyb, kzrbu, predictands, grid)
    else
        flwrbd = compute_vertical_wind(ixr, jyb, kzrbd, predictands, grid)
        flwrbu = compute_vertical_wind(ixr, jyb, kzrbu, predictands, grid)
    end

    if zrfu < ztildetfc[ixr, jyf, k0 - 1]
        flwrfd = 0.0
        flwrfu = 0.0
    elseif zrfd < ztildetfc[ixr, jyf, k0 - 1]
        flwrfd = 0.0
        flwrfu = compute_vertical_wind(ixr, jyf, kzrfu, predictands, grid)
    else
        flwrfd = compute_vertical_wind(ixr, jyf, kzrfd, predictands, grid)
        flwrfu = compute_vertical_wind(ixr, jyf, kzrfu, predictands, grid)
    end

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::DUDX,
)
    (; sizex, sizey) = namelists.domain
    (; u) = predictands
    (; io, jo, j0, k1) = domain
    (; lx, ly, dx, dy, dz, x, y, ztfc, met) = grid

    if sizex == 1
        flw = 0.0
        return flw
    else
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        if ixl - 1 < 1
            error("Error in interpolate_mean_flow (DUDX): ixl - 1 = ", ixl - 1, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error("Error in interpolate_mean_flow (DUDX): ixr = ", ixr, " > nxx = ", nxx)
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
            error("Error in interpolate_mean_flow (DUDX): jyf = ", jyf, " > nyy = ", nyy)
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

    flwlbd =
        (u[ixl, jyb, kzlbd] - u[ixl - 1, jyb, kzlbd]) / dx +
        met[ixl, jyb, kzlbd, 1, 3] *
        0.25 *
        (
            u[ixl, jyb, kzlbd + 1] + u[ixl - 1, jyb, kzlbd + 1] -
            u[ixl, jyb, kzlbd - 1] - u[ixl - 1, jyb, kzlbd - 1]
        ) / dz
    flwlbu =
        (u[ixl, jyb, kzlbu] - u[ixl - 1, jyb, kzlbu]) / dx +
        met[ixl, jyb, kzlbu, 1, 3] *
        0.25 *
        (
            u[ixl, jyb, kzlbu + 1] + u[ixl - 1, jyb, kzlbu + 1] -
            u[ixl, jyb, kzlbu - 1] - u[ixl - 1, jyb, kzlbu - 1]
        ) / dz

    flwlfd =
        (u[ixl, jyf, kzlfd] - u[ixl - 1, jyf, kzlfd]) / dx +
        met[ixl, jyf, kzlfd, 1, 3] *
        0.25 *
        (
            u[ixl, jyf, kzlfd + 1] + u[ixl - 1, jyf, kzlfd + 1] -
            u[ixl, jyf, kzlfd - 1] - u[ixl - 1, jyf, kzlfd - 1]
        ) / dz
    flwlfu =
        (u[ixl, jyf, kzlfu] - u[ixl - 1, jyf, kzlfu]) / dx +
        met[ixl, jyf, kzlfu, 1, 3] *
        0.25 *
        (
            u[ixl, jyf, kzlfu + 1] + u[ixl - 1, jyf, kzlfu + 1] -
            u[ixl, jyf, kzlfu - 1] - u[ixl - 1, jyf, kzlfu - 1]
        ) / dz

    flwrbd =
        (u[ixr, jyb, kzrbd] - u[ixr - 1, jyb, kzrbd]) / dx +
        met[ixr, jyb, kzrbd, 1, 3] *
        0.25 *
        (
            u[ixr, jyb, kzrbd + 1] + u[ixr - 1, jyb, kzrbd + 1] -
            u[ixr, jyb, kzrbd - 1] - u[ixr - 1, jyb, kzrbd - 1]
        ) / dz
    flwrbu =
        (u[ixr, jyb, kzrbu] - u[ixr - 1, jyb, kzrbu]) / dx +
        met[ixr, jyb, kzrbu, 1, 3] *
        0.25 *
        (
            u[ixr, jyb, kzrbu + 1] + u[ixr - 1, jyb, kzrbu + 1] -
            u[ixr, jyb, kzrbu - 1] - u[ixr - 1, jyb, kzrbu - 1]
        ) / dz

    flwrfd =
        (u[ixr, jyf, kzrfd] - u[ixr - 1, jyf, kzrfd]) / dx +
        met[ixr, jyf, kzrfd, 1, 3] *
        0.25 *
        (
            u[ixr, jyf, kzrfd + 1] + u[ixr - 1, jyf, kzrfd + 1] -
            u[ixr, jyf, kzrfd - 1] - u[ixr - 1, jyf, kzrfd - 1]
        ) / dz
    flwrfu =
        (u[ixr, jyf, kzrfu] - u[ixr - 1, jyf, kzrfu]) / dx +
        met[ixr, jyf, kzrfu, 1, 3] *
        0.25 *
        (
            u[ixr, jyf, kzrfu + 1] + u[ixr - 1, jyf, kzrfu + 1] -
            u[ixr, jyf, kzrfu - 1] - u[ixr - 1, jyf, kzrfu - 1]
        ) / dz

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::DUDY,
)
    (; sizex, sizey) = namelists.domain
    (; u) = predictands
    (; io, jo, i0, k1) = domain
    (; lx, ly, dx, dy, dz, x, y, ztfc, met) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DUDY): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error("Error in interpolate_mean_flow (DUDY): ixr = ", ixr, " > nxx = ", nxx)
        end
    end
    xr = x[ixr + io] + dx / 2
    xl = x[ixl + io] + dx / 2

    # Locate the closest points in meridional direction.
    if sizey == 1
        flw = 0.0
        return flw
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - jo
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

    flwlbd =
        (u[ixl, jyb + 1, kzlbd] - u[ixl, jyb, kzlbd]) / dy +
        0.25 *
        (
            met[ixl, jyb, kzlbd, 2, 3] +
            met[ixl + 1, jyb, kzlbd, 2, 3] +
            met[ixl, jyb + 1, kzlbd, 2, 3] +
            met[ixl + 1, jyb + 1, kzlbd, 2, 3]
        ) *
        0.25 *
        (
            u[ixl, jyb, kzlbd + 1] + u[ixl, jyb + 1, kzlbd + 1] -
            u[ixl, jyb, kzlbd - 1] - u[ixl, jyb + 1, kzlbd - 1]
        ) / dz
    flwlbu =
        (u[ixl, jyb + 1, kzlbu] - u[ixl, jyb, kzlbu]) / dy +
        0.25 *
        (
            met[ixl, jyb, kzlbu, 2, 3] +
            met[ixl + 1, jyb, kzlbu, 2, 3] +
            met[ixl, jyb + 1, kzlbu, 2, 3] +
            met[ixl + 1, jyb + 1, kzlbu, 2, 3]
        ) *
        0.25 *
        (
            u[ixl, jyb, kzlbu + 1] + u[ixl, jyb + 1, kzlbu + 1] -
            u[ixl, jyb, kzlbu - 1] - u[ixl, jyb + 1, kzlbu - 1]
        ) / dz

    flwlfd =
        (u[ixl, jyf + 1, kzlfd] - u[ixl, jyf, kzlfd]) / dy +
        0.25 *
        (
            met[ixl, jyf, kzlfd, 2, 3] +
            met[ixl + 1, jyf, kzlfd, 2, 3] +
            met[ixl, jyf + 1, kzlfd, 2, 3] +
            met[ixl + 1, jyf + 1, kzlfd, 2, 3]
        ) *
        0.25 *
        (
            u[ixl, jyf, kzlfd + 1] + u[ixl, jyf + 1, kzlfd + 1] -
            u[ixl, jyf, kzlfd - 1] - u[ixl, jyf + 1, kzlfd - 1]
        ) / dz
    flwlfu =
        (u[ixl, jyf + 1, kzlfu] - u[ixl, jyf, kzlfu]) / dy +
        0.25 *
        (
            met[ixl, jyf, kzlfu, 2, 3] +
            met[ixl + 1, jyf, kzlfu, 2, 3] +
            met[ixl, jyf + 1, kzlfu, 2, 3] +
            met[ixl + 1, jyf + 1, kzlfu, 2, 3]
        ) *
        0.25 *
        (
            u[ixl, jyf, kzlfu + 1] + u[ixl, jyf + 1, kzlfu + 1] -
            u[ixl, jyf, kzlfu - 1] - u[ixl, jyf + 1, kzlfu - 1]
        ) / dz

    flwrbd =
        (u[ixr, jyb + 1, kzrbd] - u[ixr, jyb, kzrbd]) / dy +
        0.25 *
        (
            met[ixr, jyb, kzrbd, 2, 3] +
            met[ixr + 1, jyb, kzrbd, 2, 3] +
            met[ixr, jyb + 1, kzrbd, 2, 3] +
            met[ixr + 1, jyb + 1, kzrbd, 2, 3]
        ) *
        0.25 *
        (
            u[ixr, jyb, kzrbd + 1] + u[ixr, jyb + 1, kzrbd + 1] -
            u[ixr, jyb, kzrbd - 1] - u[ixr, jyb + 1, kzrbd - 1]
        ) / dz
    flwrbu =
        (u[ixr, jyb + 1, kzrbu] - u[ixr, jyb, kzrbu]) / dy +
        0.25 *
        (
            met[ixr, jyb, kzrbu, 2, 3] +
            met[ixr + 1, jyb, kzrbu, 2, 3] +
            met[ixr, jyb + 1, kzrbu, 2, 3] +
            met[ixr + 1, jyb + 1, kzrbu, 2, 3]
        ) *
        0.25 *
        (
            u[ixr, jyb, kzrbu + 1] + u[ixr, jyb + 1, kzrbu + 1] -
            u[ixr, jyb, kzrbu - 1] - u[ixr, jyb + 1, kzrbu - 1]
        ) / dz

    flwrfd =
        (u[ixr, jyf + 1, kzrfd] - u[ixr, jyf, kzrfd]) / dy +
        0.25 *
        (
            met[ixr, jyf, kzrfd, 2, 3] +
            met[ixr + 1, jyf, kzrfd, 2, 3] +
            met[ixr, jyf + 1, kzrfd, 2, 3] +
            met[ixr + 1, jyf + 1, kzrfd, 2, 3]
        ) *
        0.25 *
        (
            u[ixr, jyf, kzrfd + 1] + u[ixr, jyf + 1, kzrfd + 1] -
            u[ixr, jyf, kzrfd - 1] - u[ixr, jyf + 1, kzrfd - 1]
        ) / dz
    flwrfu =
        (u[ixr, jyf + 1, kzrfu] - u[ixr, jyf, kzrfu]) / dy +
        0.25 *
        (
            met[ixr, jyf, kzrfu, 2, 3] +
            met[ixr + 1, jyf, kzrfu, 2, 3] +
            met[ixr, jyf + 1, kzrfu, 2, 3] +
            met[ixr + 1, jyf + 1, kzrfu, 2, 3]
        ) *
        0.25 *
        (
            u[ixr, jyf, kzrfu + 1] + u[ixr, jyf + 1, kzrfu + 1] -
            u[ixr, jyf, kzrfu - 1] - u[ixr, jyf + 1, kzrfu - 1]
        ) / dz

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::DUDZ,
)
    (; sizex, sizey) = namelists.domain
    (; u) = predictands
    (; io, jo, i0, j0, k0, k1) = domain
    (; lx, ly, lz, dx, dy, dz, jac, x, y, ztildetfc) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - io
        if ixl < 1
            error("Error in interpolate_mean_flow (DUDZ): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error("Error in interpolate_mean_flow (DUDZ): ixr = ", ixr, " > nxx = ", nxx)
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
            error("Error in interpolate_mean_flow (DUDZ): jyf = ", jyf, " > nyy = ", nyy)
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

    if zlbu < ztildetfc[ixl, jyb, k0]
        flwlbd = 0.0
        flwlbu = 0.0
    elseif zlbd < ztildetfc[ixl, jyb, k0]
        flwlbd = 0.0
        flwlbu =
            (u[ixl, jyb, kzlbu + 1] - u[ixl, jyb, kzlbu]) / dz / (
                jac[ixl, jyb, kzlbu] * jac[ixl, jyb, kzlbu + 1] /
                (jac[ixl, jyb, kzlbu] + jac[ixl, jyb, kzlbu + 1]) +
                jac[ixl + 1, jyb, kzlbu] * jac[ixl + 1, jyb, kzlbu + 1] /
                (jac[ixl + 1, jyb, kzlbu] + jac[ixl + 1, jyb, kzlbu + 1])
            )
    else
        if zlbu < lz[2]
            flwlbd =
                (u[ixl, jyb, kzlbd + 1] - u[ixl, jyb, kzlbd]) / dz / (
                    jac[ixl, jyb, kzlbd] * jac[ixl, jyb, kzlbd + 1] /
                    (jac[ixl, jyb, kzlbd] + jac[ixl, jyb, kzlbd + 1]) +
                    jac[ixl + 1, jyb, kzlbd] * jac[ixl + 1, jyb, kzlbd + 1] /
                    (jac[ixl + 1, jyb, kzlbd] + jac[ixl + 1, jyb, kzlbd + 1])
                )
            flwlbu =
                (u[ixl, jyb, kzlbu + 1] - u[ixl, jyb, kzlbu]) / dz / (
                    jac[ixl, jyb, kzlbu] * jac[ixl, jyb, kzlbu + 1] /
                    (jac[ixl, jyb, kzlbu] + jac[ixl, jyb, kzlbu + 1]) +
                    jac[ixl + 1, jyb, kzlbu] * jac[ixl + 1, jyb, kzlbu + 1] /
                    (jac[ixl + 1, jyb, kzlbu] + jac[ixl + 1, jyb, kzlbu + 1])
                )
        elseif zlbd < lz[2]
            flwlbd =
                (u[ixl, jyb, kzlbd + 1] - u[ixl, jyb, kzlbd]) / dz / (
                    jac[ixl, jyb, kzlbd] * jac[ixl, jyb, kzlbd + 1] /
                    (jac[ixl, jyb, kzlbd] + jac[ixl, jyb, kzlbd + 1]) +
                    jac[ixl + 1, jyb, kzlbd] * jac[ixl + 1, jyb, kzlbd + 1] /
                    (jac[ixl + 1, jyb, kzlbd] + jac[ixl + 1, jyb, kzlbd + 1])
                )
            flwlbu = 0.0
        else
            flwlbd = 0.0
            flwlbu = 0.0
        end
    end

    if zlfu < ztildetfc[ixl, jyf, k0]
        flwlfd = 0.0
        flwlfu = 0.0
    elseif zlfd < ztildetfc[ixl, jyf, k0]
        flwlfd = 0.0
        flwlfu =
            (u[ixl, jyf, kzlfu + 1] - u[ixl, jyf, kzlfu]) / dz / (
                jac[ixl, jyf, kzlfu] * jac[ixl, jyf, kzlfu + 1] /
                (jac[ixl, jyf, kzlfu] + jac[ixl, jyf, kzlfu + 1]) +
                jac[ixl + 1, jyf, kzlfu] * jac[ixl + 1, jyf, kzlfu + 1] /
                (jac[ixl + 1, jyf, kzlfu] + jac[ixl + 1, jyf, kzlfu + 1])
            )
    else
        if zlfu < lz[2]
            flwlfd =
                (u[ixl, jyf, kzlfd + 1] - u[ixl, jyf, kzlfd]) / dz / (
                    jac[ixl, jyf, kzlfd] * jac[ixl, jyf, kzlfd + 1] /
                    (jac[ixl, jyf, kzlfd] + jac[ixl, jyf, kzlfd + 1]) +
                    jac[ixl + 1, jyf, kzlfd] * jac[ixl + 1, jyf, kzlfd + 1] /
                    (jac[ixl + 1, jyf, kzlfd] + jac[ixl + 1, jyf, kzlfd + 1])
                )
            flwlfu =
                (u[ixl, jyf, kzlfu + 1] - u[ixl, jyf, kzlfu]) / dz / (
                    jac[ixl, jyf, kzlfu] * jac[ixl, jyf, kzlfu + 1] /
                    (jac[ixl, jyf, kzlfu] + jac[ixl, jyf, kzlfu + 1]) +
                    jac[ixl + 1, jyf, kzlfu] * jac[ixl + 1, jyf, kzlfu + 1] /
                    (jac[ixl + 1, jyf, kzlfu] + jac[ixl + 1, jyf, kzlfu + 1])
                )
        elseif zlfd < lz[2]
            flwlfd =
                (u[ixl, jyf, kzlfd + 1] - u[ixl, jyf, kzlfd]) / dz / (
                    jac[ixl, jyf, kzlfd] * jac[ixl, jyf, kzlfd + 1] /
                    (jac[ixl, jyf, kzlfd] + jac[ixl, jyf, kzlfd + 1]) +
                    jac[ixl + 1, jyf, kzlfd] * jac[ixl + 1, jyf, kzlfd + 1] /
                    (jac[ixl + 1, jyf, kzlfd] + jac[ixl + 1, jyf, kzlfd + 1])
                )
            flwlfu = 0.0
        else
            flwlfd = 0.0
            flwlfu = 0.0
        end
    end

    if zrbu < ztildetfc[ixr, jyb, k0]
        flwrbd = 0.0
        flwrbu = 0.0
    elseif zrbd < ztildetfc[ixr, jyb, k0]
        flwrbd = 0.0
        flwrbu =
            (u[ixr, jyb, kzrbu + 1] - u[ixr, jyb, kzrbu]) / dz / (
                jac[ixr, jyb, kzrbu] * jac[ixr, jyb, kzrbu + 1] /
                (jac[ixr, jyb, kzrbu] + jac[ixr, jyb, kzrbu + 1]) +
                jac[ixr + 1, jyb, kzrbu] * jac[ixr + 1, jyb, kzrbu + 1] /
                (jac[ixr + 1, jyb, kzrbu] + jac[ixr + 1, jyb, kzrbu + 1])
            )
    else
        if zrbu < lz[2]
            flwrbd =
                (u[ixr, jyb, kzrbd + 1] - u[ixr, jyb, kzrbd]) / dz / (
                    jac[ixr, jyb, kzrbd] * jac[ixr, jyb, kzrbd + 1] /
                    (jac[ixr, jyb, kzrbd] + jac[ixr, jyb, kzrbd + 1]) +
                    jac[ixr + 1, jyb, kzrbd] * jac[ixr + 1, jyb, kzrbd + 1] /
                    (jac[ixr + 1, jyb, kzrbd] + jac[ixr + 1, jyb, kzrbd + 1])
                )
            flwrbu =
                (u[ixr, jyb, kzrbu + 1] - u[ixr, jyb, kzrbu]) / dz / (
                    jac[ixr, jyb, kzrbu] * jac[ixr, jyb, kzrbu + 1] /
                    (jac[ixr, jyb, kzrbu] + jac[ixr, jyb, kzrbu + 1]) +
                    jac[ixr + 1, jyb, kzrbu] * jac[ixr + 1, jyb, kzrbu + 1] /
                    (jac[ixr + 1, jyb, kzrbu] + jac[ixr + 1, jyb, kzrbu + 1])
                )
        elseif zrbd < lz[2]
            flwrbd =
                (u[ixr, jyb, kzrbd + 1] - u[ixr, jyb, kzrbd]) / dz / (
                    jac[ixr, jyb, kzrbd] * jac[ixr, jyb, kzrbd + 1] /
                    (jac[ixr, jyb, kzrbd] + jac[ixr, jyb, kzrbd + 1]) +
                    jac[ixr + 1, jyb, kzrbd] * jac[ixr + 1, jyb, kzrbd + 1] /
                    (jac[ixr + 1, jyb, kzrbd] + jac[ixr + 1, jyb, kzrbd + 1])
                )
            flwrbu = 0.0
        else
            flwrbd = 0.0
            flwrbu = 0.0
        end
    end

    if zrfu < ztildetfc[ixr, jyf, k0]
        flwrfd = 0.0
        flwrfu = 0.0
    elseif zrfd < ztildetfc[ixr, jyf, k0]
        flwrfd = 0.0
        flwrbu =
            (u[ixr, jyf, kzrbu + 1] - u[ixr, jyf, kzrbu]) / dz / (
                jac[ixr, jyf, kzrbu] * jac[ixr, jyf, kzrbu + 1] /
                (jac[ixr, jyf, kzrbu] + jac[ixr, jyf, kzrbu + 1]) +
                jac[ixr + 1, jyf, kzrbu] * jac[ixr + 1, jyf, kzrbu + 1] /
                (jac[ixr + 1, jyf, kzrbu] + jac[ixr + 1, jyf, kzrbu + 1])
            )
    else
        if zrfu < lz[2]
            flwrfd =
                (u[ixr, jyf, kzrfd + 1] - u[ixr, jyf, kzrfd]) / dz / (
                    jac[ixr, jyf, kzrfd] * jac[ixr, jyf, kzrfd + 1] /
                    (jac[ixr, jyf, kzrfd] + jac[ixr, jyf, kzrfd + 1]) +
                    jac[ixr + 1, jyf, kzrfd] * jac[ixr + 1, jyf, kzrfd + 1] /
                    (jac[ixr + 1, jyf, kzrfd] + jac[ixr + 1, jyf, kzrfd + 1])
                )
            flwrfu =
                (u[ixr, jyf, kzrfu + 1] - u[ixr, jyf, kzrfu]) / dz / (
                    jac[ixr, jyf, kzrfu] * jac[ixr, jyf, kzrfu + 1] /
                    (jac[ixr, jyf, kzrfu] + jac[ixr, jyf, kzrfu + 1]) +
                    jac[ixr + 1, jyf, kzrfu] * jac[ixr + 1, jyf, kzrfu + 1] /
                    (jac[ixr + 1, jyf, kzrfu] + jac[ixr + 1, jyf, kzrfu + 1])
                )
        elseif zrfd < lz[2]
            flwrfd =
                (u[ixr, jyf, kzrfd + 1] - u[ixr, jyf, kzrfd]) / dz / (
                    jac[ixr, jyf, kzrfd] * jac[ixr, jyf, kzrfd + 1] /
                    (jac[ixr, jyf, kzrfd] + jac[ixr, jyf, kzrfd + 1]) +
                    jac[ixr + 1, jyf, kzrfd] * jac[ixr + 1, jyf, kzrfd + 1] /
                    (jac[ixr + 1, jyf, kzrfd] + jac[ixr + 1, jyf, kzrfd + 1])
                )
            flwrfu = 0.0
        else
            flwrfd = 0.0
            flwrfu = 0.0
        end
    end

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::DVDX,
)
    (; sizex, sizey) = namelists.domain
    (; v) = predictands
    (; io, jo, i0, j0, k1) = domain
    (; lx, ly, dx, dy, dz, x, y, ztfc, met) = grid

    # Locate the closest points in zonal direction.
    if sizex == 1
        flw = 0.0
        return flw
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - io
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
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - jo
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

    flwlbd =
        (v[ixl + 1, jyb, kzlbd] - v[ixl, jyb, kzlbd]) / dx +
        0.25 *
        (
            met[ixl, jyb, kzlbd, 1, 3] +
            met[ixl + 1, jyb, kzlbd, 1, 3] +
            met[ixl, jyb + 1, kzlbd, 1, 3] +
            met[ixl + 1, jyb + 1, kzlbd, 1, 3]
        ) *
        0.25 *
        (
            v[ixl, jyb, kzlbd + 1] + v[ixl + 1, jyb, kzlbd + 1] -
            v[ixl, jyb, kzlbd - 1] - v[ixl + 1, jyb, kzlbd - 1]
        ) / dz
    flwlbu =
        (v[ixl + 1, jyb, kzlbu] - v[ixl, jyb, kzlbu]) / dx +
        0.25 *
        (
            met[ixl, jyb, kzlbu, 1, 3] +
            met[ixl + 1, jyb, kzlbu, 1, 3] +
            met[ixl, jyb + 1, kzlbu, 1, 3] +
            met[ixl + 1, jyb + 1, kzlbu, 1, 3]
        ) *
        0.25 *
        (
            v[ixl, jyb, kzlbu + 1] + v[ixl + 1, jyb, kzlbu + 1] -
            v[ixl, jyb, kzlbu - 1] - v[ixl + 1, jyb, kzlbu - 1]
        ) / dz

    flwlfd =
        (v[ixl + 1, jyf, kzlfd] - v[ixl, jyf, kzlfd]) / dx +
        0.25 *
        (
            met[ixl, jyf, kzlfd, 1, 3] +
            met[ixl + 1, jyf, kzlfd, 1, 3] +
            met[ixl, jyf + 1, kzlfd, 1, 3] +
            met[ixl + 1, jyf + 1, kzlfd, 1, 3]
        ) *
        0.25 *
        (
            v[ixl, jyf, kzlfd + 1] + v[ixl + 1, jyf, kzlfd + 1] -
            v[ixl, jyf, kzlfd - 1] - v[ixl + 1, jyf, kzlfd - 1]
        ) / dz
    flwlfu =
        (v[ixl + 1, jyf, kzlfu] - v[ixl, jyf, kzlfu]) / dx +
        0.25 *
        (
            met[ixl, jyf, kzlfu, 1, 3] +
            met[ixl + 1, jyf, kzlfu, 1, 3] +
            met[ixl, jyf + 1, kzlfu, 1, 3] +
            met[ixl + 1, jyf + 1, kzlfu, 1, 3]
        ) *
        0.25 *
        (
            v[ixl, jyf, kzlfu + 1] + v[ixl + 1, jyf, kzlfu + 1] -
            v[ixl, jyf, kzlfu - 1] - v[ixl + 1, jyf, kzlfu - 1]
        ) / dz

    flwrbd =
        (v[ixr + 1, jyb, kzrbd] - v[ixr, jyb, kzrbd]) / dx +
        0.25 *
        (
            met[ixr, jyb, kzrbd, 1, 3] +
            met[ixr + 1, jyb, kzrbd, 1, 3] +
            met[ixr, jyb + 1, kzrbd, 1, 3] +
            met[ixr + 1, jyb + 1, kzrbd, 1, 3]
        ) *
        0.25 *
        (
            v[ixr, jyb, kzrbd + 1] + v[ixr + 1, jyb, kzrbd + 1] -
            v[ixr, jyb, kzrbd - 1] - v[ixr + 1, jyb, kzrbd - 1]
        ) / dz
    flwrbu =
        (v[ixr + 1, jyb, kzrbu] - v[ixr, jyb, kzrbu]) / dx +
        0.25 *
        (
            met[ixr, jyb, kzrbu, 1, 3] +
            met[ixr + 1, jyb, kzrbu, 1, 3] +
            met[ixr, jyb + 1, kzrbu, 1, 3] +
            met[ixr + 1, jyb + 1, kzrbu, 1, 3]
        ) *
        0.25 *
        (
            v[ixr, jyb, kzrbu + 1] + v[ixr + 1, jyb, kzrbu + 1] -
            v[ixr, jyb, kzrbu - 1] - v[ixr + 1, jyb, kzrbu - 1]
        ) / dz

    flwrfd =
        (v[ixr + 1, jyf, kzrfd] - v[ixr, jyf, kzrfd]) / dx +
        0.25 *
        (
            met[ixr, jyf, kzrfd, 1, 3] +
            met[ixr + 1, jyf, kzrfd, 1, 3] +
            met[ixr, jyf + 1, kzrfd, 1, 3] +
            met[ixr + 1, jyf + 1, kzrfd, 1, 3]
        ) *
        0.25 *
        (
            v[ixr, jyf, kzrfd + 1] + v[ixr + 1, jyf, kzrfd + 1] -
            v[ixr, jyf, kzrfd - 1] - v[ixr + 1, jyf, kzrfd - 1]
        ) / dz
    flwrfu =
        (v[ixr + 1, jyf, kzrfu] - v[ixr, jyf, kzrfu]) / dx +
        0.25 *
        (
            met[ixr, jyf, kzrfu, 1, 3] +
            met[ixr + 1, jyf, kzrfu, 1, 3] +
            met[ixr, jyf + 1, kzrfu, 1, 3] +
            met[ixr + 1, jyf + 1, kzrfu, 1, 3]
        ) *
        0.25 *
        (
            v[ixr, jyf, kzrfu + 1] + v[ixr + 1, jyf, kzrfu + 1] -
            v[ixr, jyf, kzrfu - 1] - v[ixr + 1, jyf, kzrfu - 1]
        ) / dz

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::DVDY,
)
    (; sizex, sizey) = namelists.domain
    (; v) = predictands
    (; io, jo, i0, j0, k1) = domain
    (; lx, ly, dx, dy, dz, x, y, ztfc, met) = grid

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
            error("Error in interpolate_mean_flow (DVDY): ixr = ", ixr, " > nxx = ", nxx)
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        flw = 0.0
        return flw
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if jyb - 1 < 1
            error("Error in interpolate_mean_flow (DVDY): jyb - 1 = ", jyb - 1, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error("Error in interpolate_mean_flow (DVDY): jyf = ", jyf, " > nyy = ", nyy)
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

    flwlbd =
        (v[ixl, jyb, kzlbd] - v[ixl, jyb - 1, kzlbd]) / dy +
        met[ixl, jyb, kzlbd, 2, 3] *
        0.25 *
        (
            v[ixl, jyb, kzlbd + 1] + v[ixl, jyb - 1, kzlbd + 1] -
            v[ixl, jyb, kzlbd - 1] - v[ixl, jyb - 1, kzlbd - 1]
        ) / dz
    flwlbu =
        (v[ixl, jyb, kzlbu] - v[ixl, jyb - 1, kzlbu]) / dy +
        met[ixl, jyb, kzlbu, 2, 3] *
        0.25 *
        (
            v[ixl, jyb, kzlbu + 1] + v[ixl, jyb - 1, kzlbu + 1] -
            v[ixl, jyb, kzlbu - 1] - v[ixl, jyb - 1, kzlbu - 1]
        ) / dz

    flwlfd =
        (v[ixl, jyf, kzlfd] - v[ixl, jyf - 1, kzlfd]) / dy +
        met[ixl, jyf, kzlfd, 2, 3] *
        0.25 *
        (
            v[ixl, jyf, kzlfd + 1] + v[ixl, jyf - 1, kzlfd + 1] -
            v[ixl, jyf, kzlfd - 1] - v[ixl, jyf - 1, kzlfd - 1]
        ) / dz
    flwlfu =
        (v[ixl, jyf, kzlfu] - v[ixl, jyf - 1, kzlfu]) / dy +
        met[ixl, jyf, kzlfu, 2, 3] *
        0.25 *
        (
            v[ixl, jyf, kzlfu + 1] + v[ixl, jyf - 1, kzlfu + 1] -
            v[ixl, jyf, kzlfu - 1] - v[ixl, jyf - 1, kzlfu - 1]
        ) / dz

    flwrbd =
        (v[ixr, jyb, kzrbd] - v[ixr, jyb - 1, kzrbd]) / dy +
        met[ixr, jyb, kzrbd, 2, 3] *
        0.25 *
        (
            v[ixr, jyb, kzrbd + 1] + v[ixr, jyb - 1, kzrbd + 1] -
            v[ixr, jyb, kzrbd - 1] - v[ixr, jyb - 1, kzrbd - 1]
        ) / dz
    flwrbu =
        (v[ixr, jyb, kzrbu] - v[ixr, jyb - 1, kzrbu]) / dy +
        met[ixr, jyb, kzrbu, 2, 3] *
        0.25 *
        (
            v[ixr, jyb, kzrbu + 1] + v[ixr, jyb - 1, kzrbu + 1] -
            v[ixr, jyb, kzrbu - 1] - v[ixr, jyb - 1, kzrbu - 1]
        ) / dz

    flwrfd =
        (v[ixr, jyf, kzrfd] - v[ixr, jyf - 1, kzrfd]) / dy +
        met[ixr, jyf, kzrfd, 2, 3] *
        0.25 *
        (
            v[ixr, jyf, kzrfd + 1] + v[ixr, jyf - 1, kzrfd + 1] -
            v[ixr, jyf, kzrfd - 1] - v[ixr, jyf - 1, kzrfd - 1]
        ) / dz
    flwrfu =
        (v[ixr, jyf, kzrfu] - v[ixr, jyf - 1, kzrfu]) / dy +
        met[ixr, jyf, kzrfu, 2, 3] *
        0.25 *
        (
            v[ixr, jyf, kzrfu + 1] + v[ixr, jyf - 1, kzrfu + 1] -
            v[ixr, jyf, kzrfu - 1] - v[ixr, jyf - 1, kzrfu - 1]
        ) / dz

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end

function interpolate_mean_flow(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    predictants::Predictands,
    flwtype::DVDZ,
)
    (; sizex, sizey) = namelists.domain
    (; v) = predictands
    (; io, jo, i0, j0, k0, k1) = domain
    (; lx, ly, lz, dx, dy, dz, jac, x, y, ztildetfc) = grid

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
            error("Error in interpolate_mean_flow (DVDZ): ixr = ", ixr, " > nxx = ", nxx)
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    # Locate the closest points in meridional direction.
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1]) / dy) + j0 - jo
        if jyb < 1
            error("Error in interpolate_mean_flow: jyb = ", jyb, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error("Error in interpolate_mean_flow: jyf = ", jyf, " > nyy = ", nyy)
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

    if zlbu < ztildetfc[ixl, jyb, k0]
        flwlbd = 0.0
        flwlbu = 0.0
    elseif zlbd < ztildetfc[ixl, jyb, k0]
        flwlbd = 0.0
        flwlbu =
            (v[ixl, jyb, kzlbu + 1] - v[ixl, jyb, kzlbu]) / dz / (
                jac[ixl, jyb, kzlbu] * jac[ixl, jyb, kzlbu + 1] /
                (jac[ixl, jyb, kzlbu] + jac[ixl, jyb, kzlbu + 1]) +
                jac[ixl, jyb + 1, kzlbu] * jac[ixl, jyb + 1, kzlbu + 1] /
                (jac[ixl, jyb + 1, kzlbu] + jac[ixl, jyb + 1, kzlbu + 1])
            )
    else
        if zlbu < lz[2]
            flwlbd =
                (v[ixl, jyb, kzlbd + 1] - v[ixl, jyb, kzlbd]) / dz / (
                    jac[ixl, jyb, kzlbd] * jac[ixl, jyb, kzlbd + 1] /
                    (jac[ixl, jyb, kzlbd] + jac[ixl, jyb, kzlbd + 1]) +
                    jac[ixl, jyb + 1, kzlbd] * jac[ixl, jyb + 1, kzlbd + 1] /
                    (jac[ixl, jyb + 1, kzlbd] + jac[ixl, jyb + 1, kzlbd + 1])
                )
            flwlbu =
                (v[ixl, jyb, kzlbu + 1] - v[ixl, jyb, kzlbu]) / dz / (
                    jac[ixl, jyb, kzlbu] * jac[ixl, jyb, kzlbu + 1] /
                    (jac[ixl, jyb, kzlbu] + jac[ixl, jyb, kzlbu + 1]) +
                    jac[ixl, jyb + 1, kzlbu] * jac[ixl, jyb + 1, kzlbu + 1] /
                    (jac[ixl, jyb + 1, kzlbu] + jac[ixl, jyb + 1, kzlbu + 1])
                )
        elseif zlbd < lz[2]
            flwlbd =
                (v[ixl, jyb, kzlbd + 1] - v[ixl, jyb, kzlbd]) / dz / (
                    jac[ixl, jyb, kzlbd] * jac[ixl, jyb, kzlbd + 1] /
                    (jac[ixl, jyb, kzlbd] + jac[ixl, jyb, kzlbd + 1]) +
                    jac[ixl, jyb + 1, kzlbd] * jac[ixl, jyb + 1, kzlbd + 1] /
                    (jac[ixl, jyb + 1, kzlbd] + jac[ixl, jyb + 1, kzlbd + 1])
                )
            flwlbu = 0.0
        else
            flwlbd = 0.0
            flwlbu = 0.0
        end
    end

    if zlfu < ztildetfc[ixl, jyf, k0]
        flwlfd = 0.0
        flwlfu = 0.0
    elseif zlfd < ztildetfc[ixl, jyf, k0]
        flwlfd = 0.0
        flwlfu =
            (v[ixl, jyf, kzlfu + 1] - v[ixl, jyf, kzlfu]) / dz / (
                jac[ixl, jyf, kzlfu] * jac[ixl, jyf, kzlfu + 1] /
                (jac[ixl, jyf, kzlfu] + jac[ixl, jyf, kzlfu + 1]) +
                jac[ixl, jyf + 1, kzlfu] * jac[ixl, jyf + 1, kzlfu + 1] /
                (jac[ixl, jyf + 1, kzlfu] + jac[ixl, jyf + 1, kzlfu + 1])
            )
    else
        if zlfu < lz[2]
            flwlfd =
                (v[ixl, jyf, kzlfd + 1] - v[ixl, jyf, kzlfd]) / dz / (
                    jac[ixl, jyf, kzlfd] * jac[ixl, jyf, kzlfd + 1] /
                    (jac[ixl, jyf, kzlfd] + jac[ixl, jyf, kzlfd + 1]) +
                    jac[ixl, jyf + 1, kzlfd] * jac[ixl, jyf + 1, kzlfd + 1] /
                    (jac[ixl, jyf + 1, kzlfd] + jac[ixl, jyf + 1, kzlfd + 1])
                )
            flwlfu =
                (v[ixl, jyf, kzlfu + 1] - v[ixl, jyf, kzlfu]) / dz / (
                    jac[ixl, jyf, kzlfu] * jac[ixl, jyf, kzlfu + 1] /
                    (jac[ixl, jyf, kzlfu] + jac[ixl, jyf, kzlfu + 1]) +
                    jac[ixl, jyf + 1, kzlfu] * jac[ixl, jyf + 1, kzlfu + 1] /
                    (jac[ixl, jyf + 1, kzlfu] + jac[ixl, jyf + 1, kzlfu + 1])
                )
        elseif zlfd < lz[2]
            flwlfd =
                (v[ixl, jyf, kzlfd + 1] - v[ixl, jyf, kzlfd]) / dz / (
                    jac[ixl, jyf, kzlfd] * jac[ixl, jyf, kzlfd + 1] /
                    (jac[ixl, jyf, kzlfd] + jac[ixl, jyf, kzlfd + 1]) +
                    jac[ixl, jyf + 1, kzlfd] * jac[ixl, jyf + 1, kzlfd + 1] /
                    (jac[ixl, jyf + 1, kzlfd] + jac[ixl, jyf + 1, kzlfd + 1])
                )
            flwlfu = 0.0
        else
            flwlfd = 0.0
            flwlfu = 0.0
        end
    end

    if zrbu < ztildetfc[ixr, jyb, k0]
        flwrbd = 0.0
        flwrbu = 0.0
    elseif zrbd < ztildetfc[ixr, jyb, k0]
        flwrbd = 0.0
        flwrbu =
            (v[ixr, jyb, kzrbu + 1] - v[ixr, jyb, kzrbu]) / dz / (
                jac[ixr, jyb, kzrbu] * jac[ixr, jyb, kzrbu + 1] /
                (jac[ixr, jyb, kzrbu] + jac[ixr, jyb, kzrbu + 1]) +
                jac[ixr, jyb + 1, kzrbu] * jac[ixr, jyb + 1, kzrbu + 1] /
                (jac[ixr, jyb + 1, kzrbu] + jac[ixr, jyb + 1, kzrbu + 1])
            )
    else
        if zrbu < lz[2]
            flwrbd =
                (v[ixr, jyb, kzrbd + 1] - v[ixr, jyb, kzrbd]) / dz / (
                    jac[ixr, jyb, kzrbd] * jac[ixr, jyb, kzrbd + 1] /
                    (jac[ixr, jyb, kzrbd] + jac[ixr, jyb, kzrbd + 1]) +
                    jac[ixr, jyb + 1, kzrbd] * jac[ixr, jyb + 1, kzrbd + 1] /
                    (jac[ixr, jyb + 1, kzrbd] + jac[ixr, jyb + 1, kzrbd + 1])
                )
            flwrbu =
                (v[ixr, jyb, kzrbu + 1] - v[ixr, jyb, kzrbu]) / dz / (
                    jac[ixr, jyb, kzrbu] * jac[ixr, jyb, kzrbu + 1] /
                    (jac[ixr, jyb, kzrbu] + jac[ixr, jyb, kzrbu + 1]) +
                    jac[ixr, jyb + 1, kzrbu] * jac[ixr, jyb + 1, kzrbu + 1] /
                    (jac[ixr, jyb + 1, kzrbu] + jac[ixr, jyb + 1, kzrbu + 1])
                )
        elseif zrbd < lz[2]
            flwrbd =
                (v[ixr, jyb, kzrbd + 1] - v[ixr, jyb, kzrbd]) / dz / (
                    jac[ixr, jyb, kzrbd] * jac[ixr, jyb, kzrbd + 1] /
                    (jac[ixr, jyb, kzrbd] + jac[ixr, jyb, kzrbd + 1]) +
                    jac[ixr, jyb + 1, kzrbd] * jac[ixr, jyb + 1, kzrbd + 1] /
                    (jac[ixr, jyb + 1, kzrbd] + jac[ixr, jyb + 1, kzrbd + 1])
                )
            flwrbu = 0.0
        else
            flwrbd = 0.0
            flwrbu = 0.0
        end
    end

    if zrfu < ztildetfc[ixr, jyf, k0]
        flwrfd = 0.0
        flwrfu = 0.0
    elseif zrfd < ztildetfc[ixr, jyf, k0]
        flwrfd = 0.0
        flwrfu =
            (v[ixr, jyf, kzrfu + 1] - v[ixr, jyf, kzrfu]) / dz / (
                jac[ixr, jyf, kzrfu] * jac[ixr, jyf, kzrfu + 1] /
                (jac[ixr, jyf, kzrfu] + jac[ixr, jyf, kzrfu + 1]) +
                jac[ixr, jyf + 1, kzrfu] * jac[ixr, jyf + 1, kzrfu + 1] /
                (jac[ixr, jyf + 1, kzrfu] + jac[ixr, jyf + 1, kzrfu + 1])
            )
    else
        if zrfu < lz[2]
            flwrfd =
                (v[ixr, jyf, kzrfd + 1] - v[ixr, jyf, kzrfd]) / dz / (
                    jac[ixr, jyf, kzrfd] * jac[ixr, jyf, kzrfd + 1] /
                    (jac[ixr, jyf, kzrfd] + jac[ixr, jyf, kzrfd + 1]) +
                    jac[ixr, jyf + 1, kzrfd] * jac[ixr, jyf + 1, kzrfd + 1] /
                    (jac[ixr, jyf + 1, kzrfd] + jac[ixr, jyf + 1, kzrfd + 1])
                )
            flwrfu =
                (v[ixr, jyf, kzrfu + 1] - v[ixr, jyf, kzrfu]) / dz / (
                    jac[ixr, jyf, kzrfu] * jac[ixr, jyf, kzrfu + 1] /
                    (jac[ixr, jyf, kzrfu] + jac[ixr, jyf, kzrfu + 1]) +
                    jac[ixr, jyf + 1, kzrfu] * jac[ixr, jyf + 1, kzrfu + 1] /
                    (jac[ixr, jyf + 1, kzrfu] + jac[ixr, jyf + 1, kzrfu + 1])
                )
        elseif zrfd < lz[2]
            flwrfd =
                (v[ixr, jyf, kzrfd + 1] - v[ixr, jyf, kzrfd]) / dz / (
                    jac[ixr, jyf, kzrfd] * jac[ixr, jyf, kzrfd + 1] /
                    (jac[ixr, jyf, kzrfd] + jac[ixr, jyf, kzrfd + 1]) +
                    jac[ixr, jyf + 1, kzrfd] * jac[ixr, jyf + 1, kzrfd + 1] /
                    (jac[ixr, jyf + 1, kzrfd] + jac[ixr, jyf + 1, kzrfd + 1])
                )
            flwrfu = 0.0
        else
            flwrfd = 0.0
            flwrfu = 0.0
        end
    end

    # Interpolate.
    flw = interpolate(
        namelists,
        flwlbd,
        flwlbu,
        flwlfd,
        flwlfu,
        flwrbd,
        flwrbu,
        flwrfd,
        flwrfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zlc,
        xl,
        xr,
        xlc,
        yf,
        yb,
        ylc,
    )

    return flw
end
