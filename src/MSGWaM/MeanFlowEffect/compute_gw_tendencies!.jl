function compute_gw_tendencies!(state)
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; zmin_wkb_dim) = state.namelists.wkb
    (; tref, lref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; lz, dx, dy, dz, ztfc, jac, met) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (; integrals, tendencies) = state.wkb

    # Set the Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(GWTendencies)
        getfield(tendencies, field) .= 0.0
    end

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if ztfc[ix, jy, kz] < lz[1] + zmin_wkb_dim / lref
            continue
        end

        rhotot = rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]

        # Compute the drag on the zonal wind.

        tendencies.dudt[ix, jy, kz] =
            -rhotot / rhostrattfc[ix, jy, kz] / jac[ix, jy, kz] *
            (integrals.uw[ix, jy, kz + 1] - integrals.uw[ix, jy, kz - 1]) /
            (2.0 * dz)

        if sizex > 1
            tendencies.dudt[ix, jy, kz] -=
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
            tendencies.dudt[ix, jy, kz] -=
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

        tendencies.dudt[ix, jy, kz] += rhotot * integrals.etx[ix, jy, kz]

        # Compute the drag on the meridional wind.

        tendencies.dvdt[ix, jy, kz] =
            -rhotot / rhostrattfc[ix, jy, kz] / jac[ix, jy, kz] *
            (integrals.vw[ix, jy, kz + 1] - integrals.vw[ix, jy, kz - 1]) /
            (2.0 * dz)

        if sizex > 1
            tendencies.dvdt[ix, jy, kz] -=
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
            tendencies.dvdt[ix, jy, kz] -=
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

        tendencies.dvdt[ix, jy, kz] += rhotot * integrals.ety[ix, jy, kz]

        # Compute the heating.

        if fc != 0.0 && (sizex > 1 || sizey > 1)
            if sizex > 1
                integrals.dthetadt[ix, jy, kz] +=
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
                integrals.dthetadt[ix, jy, kz] +=
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
