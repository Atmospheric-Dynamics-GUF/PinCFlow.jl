function calc_gw_tendencies!(state)

    # wave impact on horizontal momentum and wave-induced heating
    # integrals.dudt, dvdt, dthetadt

    (; rho) = state.variables.predictands
    (; rhostrattfc) = state.atmosphere
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dx, dy, dz) = state.grid
    (; sizex, sizey) = state.namelists.domain
    (; integrals) = state.wkb

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        rhotot = rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]

        # forcing in x direction

        integrals.dudt[ix, jy, kz] =
            -rhotot / rhostrattfc[ix, jy, kz] / jac[ix, jy, kz] *
            (integrals.uw[ix, jy, kz + 1] - integrals.uw[ix, jy, kz - 1]) /
            (2.0 * dz)

        if sizex > 1
            integrals.dudt[ix, jy, kz] =
                integrals.dudt[ix, jy, kz] -
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.uu[ix + 1, jy, kz] -
                        integrals.uu[ix - 1, jy, kz]
                    ) / (2.0 * dx) +
                    met[ix, jy, kz, 1, 3] * (
                        integrals.uu[ix, jy, kz + 1] -
                        integrals.uu[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        if sizey > 1
            integrals.dudt[ix, jy, kz] =
                integrals.dudt[ix, jy, kz] -
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.uv[ix, jy + 1, kz] -
                        integrals.uv[ix, jy - 1, kz]
                    ) / (2.0 * dy) +
                    met[ix, jy, kz, 2, 3] * (
                        integrals.uv[ix, jy, kz + 1] -
                        integrals.uv[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        integrals.dudt[ix, jy, kz] =
            integrals.dudt[ix, jy, kz] + rhotot * integrals.etx[ix, jy, kz]

        # forcing in y direction

        rhotot = rho[ix, jy, kz]
        rhotot = rhotot + rhostrattfc[ix, jy, kz]

        integrals.dvdt[ix, jy, kz] =
            -rhotot / rhostrattfc[ix, jy, kz] / jac[ix, jy, kz] *
            (integrals.vw[ix, jy, kz + 1] - integrals.vw[ix, jy, kz - 1]) /
            (2.0 * dz)

        if sizex > 1
            integrals.dvdt[ix, jy, kz] =
                integrals.dvdt[ix, jy, kz] -
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.uv[ix + 1, jy, kz] -
                        integrals.uv[ix - 1, jy, kz]
                    ) / (2.0 * dx) +
                    met[ix, jy, kz, 1, 3] * (
                        integrals.uv[ix, jy, kz + 1] -
                        integrals.uv[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        if sizey > 1
            integrals.dvdt[ix, jy, kz] =
                integrals.dvdt[ix, jy, kz] -
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.vv[ix, jy + 1, kz] -
                        integrals.vv[ix, jy - 1, kz]
                    ) / (2.0 * dy) +
                    met[ix, jy, kz, 2, 3] * (
                        integrals.vv[ix, jy, kz + 1] -
                        integrals.vv[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        integrals.dvdt[ix, jy, kz] =
            integrals.dvdt[ix, jy, kz] + rhotot * integrals.ety[ix, jy, kz]
    end

    if (f_cor_nd != 0.0) && (sizex > 1 || sizey > 1)
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            rhotot = var.rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]

            if sizex > 1
                integrals.dthetadt[ix, jy, kz] =
                    rhotot * (
                        (
                            integrals.utheta[ix + 1, jy, kz] -
                            integrals.utheta[ix - 1, jy, kz]
                        ) / (2.0 * dx) +
                        met[ix, jy, kz, 1, 3] * (
                            integrals.utheta[ix, jy, kz + 1] -
                            integrals.utheta[ix, jy, kz - 1]
                        ) / (2.0 * dz)
                    )
            end

            if sizey > 1
                integrals.dthetadt[ix, jy, kz] =
                    integrals.dthetadt[ix, jy, kz] +
                    rhotot * (
                        (
                            integrals.vtheta[ix, jy + 1, kz] -
                            integrals.vtheta[ix, jy - 1, kz]
                        ) / (2.0 * dy) +
                        met[ix, jy, kz, 2, 3] * (
                            integrals.vtheta[ix, jy, kz + 1] -
                            integrals.vtheta[ix, jy, kz - 1]
                        ) / (2.0 * dz)
                    )
            end
        end
    end
end
