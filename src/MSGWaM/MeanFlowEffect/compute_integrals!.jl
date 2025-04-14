function compute_integrals!(state::State, mode::MultiColumn)
    (; domain, grid) = state
    (; sizex, sizey, nbx, nby, nbz) = state.namelists.domain
    (; f_coriolis_dim) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; i0, i1, j0, j1, k0, k1, io, jo) = domain
    (; lx, ly, dx, dy, dz, ztildetfc, jac) = grid
    (; rhostrattfc, thetastrattfc) = state.atmosphere
    (; nray, rays, integrals) = state.wkb

    # Set Coriolis parameter.
    f_cor_nd = f_coriolis_dim * tref

    set_integrals_to_zero!(integrals)

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

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

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            # squared Brunt-Vaisala frequency
            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + f_cor_nd^2 * mr^2) /
                sqrt(khr^2 + mr^2)

            cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            cgirz = -mr * (omir^2 - f_cor_nd^2) / (omir * (khr^2 + mr^2))

            (ixmin, ixmax, jymin, jymax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min((xr + dxr / 2), (lx[1] + (ix + io) * dx)) -
                        max((xr - dxr / 2), (lx[1] + (ix + io - 1) * dx))
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyi = (
                            min((yr + dyr / 2), (ly[1] + (jy + jo) * dy)) -
                            max((yr - dyr / 2), (ly[1] + (jy + jo - 1) * dy))
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        if sizex > 1
                            if f_cor_nd != 0.0
                                integrals.uu[ix, jy, kz] +=
                                    wadr * (
                                        kr * cgirx -
                                        (kr * cgirx + lr * cgiry) /
                                        (1.0 - (omir / f_cor_nd)^2)
                                    )
                            else
                                integrals.uu[ix, jy, kz] += wadr * kr * cgirx
                            end
                        end

                        if sizex > 1 || sizey > 1
                            integrals.uv[ix, jy, kz] += wadr * cgirx * lr
                        end

                        integrals.uw[ix, jy, kz] +=
                            wadr * kr * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        if sizey > 1
                            if f_cor_nd != 0.0
                                integrals.vv[ix, jy, kz] +=
                                    wadr * (
                                        lr * cgiry -
                                        (kr * cgirx + lr * cgiry) /
                                        (1.0 - (omir / f_cor_nd)^2)
                                    )
                            else
                                integrals.vv[ix, jy, kz] += wadr * lr * cgiry
                            end
                        end

                        integrals.vw[ix, jy, kz] +=
                            wadr * lr * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        if f_cor_nd != 0.0
                            integrals.etx[ix, jy, kz] +=
                                wadr * f_cor_nd^2 * n2r * kr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )

                            integrals.ety[ix, jy, kz] +=
                                wadr * f_cor_nd^2 * n2r * lr * mr / (
                                    rhostrattfc[ix, jy, kz] *
                                    g_ndim *
                                    omir *
                                    (khr^2 + mr^2)
                                )
                        end

                        integrals.e[ix, jy, kz] += wadr * omir
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

function compute_integrals!(
    state::State,
    mode::Union{SingleColumn, SteadyState},
)

    # steady_state or single_column
    # only calculate integrals.uw, vw, and e
    (; domain, grid) = state
    (; nbx, nby, nbz, i0, i1, j0, j1, k0, k1, io, jo) = state.domain
    (; lx, ly, dx, dy, dz, ztildetfc, jac) = state.grid
    (; sizex, sizey) = state.namelists.domain
    (; branchr) = state.namelists.wkb
    (; nray, rays, integrals) = state.wkb
    (; f_cor_nd) = state.atmosphere

    set_integrals_to_zero!(integrals)

    for kzrv in (k0 - 1):(k1 + 1),
        jyrv in (j0 - 1):(j1 + 1),
        ixrv in (i0 - 1):(i1 + 1)

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

            kr = rays.k[iray, ixrv, jyrv, kzrv]
            lr = rays.l[iray, ixrv, jyrv, kzrv]
            mr = rays.m[iray, ixrv, jyrv, kzrv]

            dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
            dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
            dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

            khr = sqrt(kr^2 + lr^2)

            # squared Brunt-Vaisala frequency
            n2r = interpolate_stratification(zr, state, N2())

            omir =
                branchr * sqrt(n2r * khr^2 + f_cor_nd^2 * mr^2) /
                sqrt(khr^2 + mr^2)

            cgirz = -mr * (omir^2 - f_cor_nd^2) / (omir * (khr^2 + mr^2))

            ixmin, ixmax, jymin, jymax =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            for ix in ixmin:ixmax
                if sizex > 1
                    dxi = (
                        min((xr + dxr / 2), (lx[1] + (ix + io) * dx)) -
                        max((xr - dxr / 2), (lx[1] + (ix + io - 1) * dx))
                    )

                    fcpspx = dkr * dxi / dx
                else
                    fcpspx = 1.0
                end

                for jy in jymin:jymax
                    if sizey > 1
                        dyr = (
                            min((yr + dyr / 2), (ly[1] + (jy + jo) * dy)) -
                            max((yr - dyr / 2), (ly[1] + (jy + jo - 1) * dy))
                        )

                        fcpspy = dlr * dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kzmin =
                        get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
                    kzmax =
                        get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)

                    for kz in kzmin:kzmax
                        dzi =
                            min((zr + dzr / 2), ztildetfc[ix, jy, kz]) -
                            max((zr - dzr / 2), ztildetfc[ix, jy, kz - 1])

                        fcpspz = dmr * dzi / jac[ix, jy, kz] / dz

                        wadr =
                            fcpspx *
                            fcpspy *
                            fcpspz *
                            rays.dens[iray, ixrv, jyrv, kzrv]

                        integrals.uw[ix, jy, kz] +=
                            wadr * kr * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        integrals.vw[ix, jy, kz] +=
                            wadr * lr * cgirz / (1.0 - (f_cor_nd / omir)^2)

                        integrals.e[ix, jy, kz] += wadr * omir
                    end
                end
            end
        end
    end
end
