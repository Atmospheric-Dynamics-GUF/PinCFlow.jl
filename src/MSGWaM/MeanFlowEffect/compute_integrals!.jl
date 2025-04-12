function compute_integrals!(state::State, mode::MultiColumn)
    (; domain, grid) = state
    (; nbx, nby, nbz, i0, i1, j0, j1, k0, k1, io, jo) = state.domain
    (; lx, ly, dx, dy, dz, ztilfetfc, jac) = state.grid
    (; sizex, sizey) = state.namelist.domain
    (; branchr) = state.namelist.WKBNamelist
    (; rays, integrals) = state.WKB
    (; f_cor_nd, thetastrattfc) = state.atmosphere

    set_integrals_to_zero!(integrals)

    for kzrv in (k0 - nbz):(k1 + nbz),
        jyrv in (j0 - nby):(j1 + nby),
        ixrv in (i0 - nbx):(i1 + nbx)

        if nray[ixrv, jyrv, kzrv] < 1
            continue
        end

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0.0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            wnrk = rays.k[iray, ixrv, jyrv, kzrv]
            wnrl = rays.l[iray, ixrv, jyrv, kzrv]
            wnrm = rays.m[iray, ixrv, jyrv, kzrv]

            dwnrk = rays.dkray[iray, ixrv, jyrv, kzrv]
            dwnrl = rays.dlray[iray, ixrv, jyrv, kzrv]
            dwnrm = rays.dmray[iray, ixrv, jyrv, kzrv]

            wnrh = sqrt(wnrk^2.0 + wnrl^2.0)

            # squared Brunt-Vaisala frequency
            nnr = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(nnr * wnrh^2.0 + f_cor_nd^2.0 * wnrm^2.0) /
                sqrt(wnrh^2.0 + wnrm^2.0)

            cgirx = wnrk * (nnr - omir^2.0) / (omir * (wnrh^2.0 + wnrm^2.0))
            cgiry = wnrl * (nnr - omir^2.0) / (omir * (wnrh^2.0 + wnrm^2.0))
            cgirz =
                -wnrm * (omir^2.0 - f_cor_nd^2.0) /
                (omir * (wnrh^2.0 + wnrm^2.0))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min((xr + dxr * 0.5), (lx[0] + (ix + io) * dx)) -
                        max((xr - dxr * 0.5), (lx[0] + (ix + io - 1) * dx))
                    )

                    fcpspx = dwnrk * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyr = (
                            min((yr + dyr * 0.5), (ly[0] + (jy + jo) * dy)) - max(
                                (yr - dyr * 0.5),
                                (ly[0] + (jy + jo - 1) * dy),
                            )
                        )

                        fcpspy = dwnrl * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin = get_next_half_level(
                        ix,
                        jy,
                        zr - dzr * 0.5,
                        domain,
                        grid,
                    )
                    kzmax = get_next_half_level(
                        ix,
                        jy,
                        zr + dzr * 0.5,
                        domain,
                        grid,
                    )

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr * 0.5), ztilfetfc[ix, jy, kz]) -
                            max((zr - dzr * 0.5), ztilfetfc[ix, jy, kz])

                        fcpspz = dwnrm * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, iyrv, izrv]

                        if sizex > 1
                            if f_cor_nd != 0.0
                                integrals.uu[ix, jy, kz] =
                                    integrals.uu[ix, jy, kz] +
                                    wadr * (
                                        wnrk * cgirx -
                                        (wnrk * cgirx + wnrl * cgiry) /
                                        (1.0 - (omir / f_cor_nd)^2)
                                    )
                            else
                                integrals.uu[ix, jy, kz] =
                                    integrals.uu[ix, jy, kz] +
                                    wadr * wnrk * cgirx
                            end
                        end

                        if sizex > 1 || sizey > 1
                            integrals.uv[ix, jy, kz] =
                                integrals.uv[ix, jy, kz] + wadr * cgirx * wnrl
                        end

                        integrals.uw[ix, jy, kz] =
                            integrals.uw[ix, jy, kz] +
                            wadr * wnrk * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        if sizey > 1
                            if f_cor_nd != 0.0
                                integrals.vv[ix, jy, kz] =
                                    integrals.vv[ix, jy, kz] +
                                    wadr * (
                                        wnrl * cgiry -
                                        (wnrk * cgirx + wnrl * cgiry) /
                                        (1.0 - (omir / f_cor_nd)^2)
                                    )
                            else
                                integrals.vv[ix, jy, kz] =
                                    integrals.vv[ix, jy, kz] +
                                    wadr * wnrl * cgiry
                            end
                        end

                        integrals.vw[ix, jy, kz] =
                            integrals.vw[ix, jy, kz] +
                            wadr * wnrl * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        if f_cor_nd != 0.0
                            integrals.etx[ix, jy, kz] =
                                integrals.etx[ix, jy, kz] +
                                wadr * f_cor_nd^2 * nnr * wnrk * wnrm / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (wnrh^2 + wnrm^2)
                                )

                            integrals.ety[ix, jy, kz] =
                                integrals.ety[ix, jy, kz] +
                                wadr * f_cor_nd^2 * nnr * wnrl * wnrm / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (wnrh^2 + wnrm^2)
                                )
                        end

                        integrals.e[ix, jy, kz] =
                            integrals.e[ix, jy, kz] + wadr * omir
                    end
                end
            end
        end
    end

    if f_cor_nd != 0.0
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            integrals.utheta[ix, jy, kz] =
                thetastrattfc[ix, jy, kz] / f_cor_nd * integrals.ety[ix, jy, kz]
            integrals.vtheta[ix, jy, kz] =
                -thetastrattfc[ix, jy, kz] / f_cor_nd *
                integrals.etx[ix, jy, kz]
        end
    end
end

function compute_integrals!(state::State, mode::AbstractWKBMode)

    # steady_state or single_column
    # only calculate integrals.uw, vw, and e
    (; domain, grid) = state
    (; nbx, nby, nbz, i0, i1, j0, j1, k0, k1, io, jo) = state.domain
    (; lx, ly, dx, dy, dz, ztilfetfc, jac) = state.grid
    (; sizex, sizey) = state.namelist.domain
    (; branchr) = state.namelist.WKBNamelist
    (; rays, integrals) = state.WKB
    (; f_cor_nd) = state.atmosphere

    set_integrals_to_zero!(integrals)

    for kzrv in (k0 - nbz):(k1 + nbz),
        jyrv in (j0 - nby):(j1 + nby),
        ixrv in (i0 - nbx):(i1 + nbx)

        if nray[ixrv, jyrv, kzrv] < 1
            continue
        end

        for iray in 1:nray[ixrv, jyrv, kzrv]
            if rays.dens[iray, ixrv, jyrv, kzrv] == 0.0
                continue
            end

            xr = rays.x[iray, ixrv, jyrv, kzrv]
            yr = rays.y[iray, ixrv, jyrv, kzrv]
            zr = rays.z[iray, ixrv, jyrv, kzrv]

            dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
            dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
            dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

            wnrk = rays.k[iray, ixrv, jyrv, kzrv]
            wnrl = rays.l[iray, ixrv, jyrv, kzrv]
            wnrm = rays.m[iray, ixrv, jyrv, kzrv]

            dwnrk = rays.dkray[iray, ixrv, jyrv, kzrv]
            dwnrl = rays.dlray[iray, ixrv, jyrv, kzrv]
            dwnrm = rays.dmray[iray, ixrv, jyrv, kzrv]

            wnrh = sqrt(wnrk^2.0 + wnrl^2.0)

            # squared Brunt-Vaisala frequency
            nnr = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(nnr * wnrh^2.0 + f_cor_nd^2.0 * wnrm^2.0) /
                sqrt(wnrh^2.0 + wnrm^2.0)

            cgirz =
                -wnrm * (omir^2.0 - f_cor_nd^2.0) /
                (omir * (wnrh^2.0 + wnrm^2.0))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min((xr + dxr * 0.5), (lx[0] + (ix + io) * dx)) -
                        max((xr - dxr * 0.5), (lx[0] + (ix + io - 1) * dx))
                    )

                    fcpspx = dwnrk * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyr = (
                            min((yr + dyr * 0.5), (ly[0] + (jy + jo) * dy)) - max(
                                (yr - dyr * 0.5),
                                (ly[0] + (jy + jo - 1) * dy),
                            )
                        )

                        fcpspy = dwnrl * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin = get_next_half_level(
                        ix,
                        jy,
                        zr - dzr * 0.5,
                        domain,
                        grid,
                    )
                    kzmax = get_next_half_level(
                        ix,
                        jy,
                        zr + dzr * 0.5,
                        domain,
                        grid,
                    )

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr * 0.5), ztilfetfc[ix, jy, kz]) -
                            max((zr - dzr * 0.5), ztilfetfc[ix, jy, kz])

                        fcpspz = dwnrm * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, iyrv, izrv]

                        integrals.uw[ix, jy, kz] =
                            integrals.uw[ix, jy, kz] +
                            wadr * wnrk * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        integrals.vw[ix, jy, kz] =
                            integrals.vw[ix, jy, kz] +
                            wadr * wnrl * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        integrals.e[ix, jy, kz] =
                            integrals.e[ix, jy, kz] + wadr * omir
                    end
                end
            end
        end
    end
end
