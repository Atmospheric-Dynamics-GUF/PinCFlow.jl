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
    (; nbz) = state.namelists.domain
    (; g_ndim) = state.constants
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; dz, jac) = state.grid
    (; pstrattfc, bvsstrattfc, thetastrattfc, rhostrattfc) = state.atmosphere
    (; rho, p) = predictands

    pstrattfc .= p

    kz0 = ko == 0 ? k0 : 1
    kz1 = ko + nzz == sizezz ? k1 : nzz

    if ko == 0
        for k in 1:nbz
            bvsstrattfc[:, :, k] .=
                g_ndim .* p[:, :, k0 - 1] ./
                (rho[:, :, k0 - 1] .+ rhostrattfc[:, :, k0 - 1]) ./
                (thetastrattfc[:, :, k0 - 1] .^ 2) ./ jac[:, :, k0 - 1] .*
                (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
        end
    end

    for kz in kz0:kz1
        bvsstrattfc[:, :, kz] .=
            g_ndim .* p[:, :, kz] ./ (rho[:, :, kz] .+ rhostrattfc[:, :, kz]) ./
            (thetastrattfc[:, :, kz] .^ 2) ./ jac[:, :, kz] .*
            (thetastrattfc[:, :, kz + 1] .- thetastrattfc[:, :, kz - 1]) ./ 2 ./
            dz
    end

    if ko + nzz == sizezz
        for k in 1:nbz
        bvsstrattfc[:, :, k1 + k] .=
            g_ndim .* p[:, :, k1 + 1] ./
            (rho[:, :, k1 + 1] .+ rhostrattfc[:, :, k1 + 1]) ./
            (thetastrattfc[:, :, k1 + 1] .^ 2) ./ jac[:, :, k1 + 1] .*
            (thetastrattfc[:, :, k1 + 1] .- thetastrattfc[:, :, k1]) ./ dz
        end
    end

    return
end
