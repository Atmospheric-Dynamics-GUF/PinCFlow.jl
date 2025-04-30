function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
)
    (; model) = state.namelists.setting
    synchronize_compressible_atmosphere!(state, predictands, model)
    return
end

function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
)
    return
end

function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::Compressible,
)
    (; g_ndim) = state.constants
    (; nxx, nyy, k0, k1) = state.domain
    (; dz, jac) = state.grid
    (; pstrattfc, bvsstrattfc, thetastrattfc, rhostrattfc) = state.atmosphere
    (; rho, p) = predictands

    pstrattfc .= p

    for jy in 1:nyy, ix in 1:nxx
        bvsstrattfc[ix, jy, k0 - 2] =
            g_ndim * p[ix, jy, k0 - 1] /
            (rho[ix, jy, k0 - 1] + rhostrattfc[ix, jy, k0 - 1]) /
            (thetastrattfc[ix, jy, k0 - 1]^2) / jac[ix, jy, k0 - 1] *
            (thetastrattfc[ix, jy, k0] - thetastrattfc[ix, jy, k0 - 1]) / dz
        bvsstrattfc[ix, jy, k0 - 1] = bvsstrattfc[ix, jy, k0 - 2]

        for kz in k0:k1
            bvsstrattfc[ix, jy, kz] =
                g_ndim * p[ix, jy, kz] /
                (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) /
                (thetastrattfc[ix, jy, kz]^2) / jac[ix, jy, kz] * (
                    thetastrattfc[ix, jy, kz + 1] -
                    thetastrattfc[ix, jy, kz - 1]
                ) / 2 / dz
        end

        bvsstrattfc[ix, jy, k1 + 1] =
            g_ndim * p[ix, jy, k1 + 1] /
            (rho[ix, jy, k1 + 1] + rhostrattfc[ix, jy, k1 + 1]) /
            (thetastrattfc[ix, jy, k1 + 1]^2) / jac[ix, jy, k1 + 1] *
            (thetastrattfc[ix, jy, k1 + 1] - thetastrattfc[ix, jy, k1]) / dz
        bvsstrattfc[ix, jy, k1 + 2] = bvsstrattfc[ix, jy, k1 + 1]
    end

    return
end
