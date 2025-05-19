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
    (; namelists, domain) = state
    (; nbz) = namelists.domain
    (; zboundaries) = namelists.setting
    (; g_ndim) = state.constants
    (; sizezz, nzz, ko, k0, k1) = domain
    (; dz, jac) = state.grid
    (; pstrattfc, bvsstrattfc, thetastrattfc, rhostrattfc) = state.atmosphere
    (; rho, p) = predictands

    pstrattfc .= p

    for kz in k0:k1
        bvsstrattfc[:, :, kz] .=
            g_ndim .* p[:, :, kz] ./ (rho[:, :, kz] .+ rhostrattfc[:, :, kz]) ./
            (thetastrattfc[:, :, kz] .^ 2) ./ jac[:, :, kz] .*
            (thetastrattfc[:, :, kz + 1] .- thetastrattfc[:, :, kz - 1]) ./ 2 ./
            dz
    end

    set_vertical_boundaries_of_field!(
        bvsstrattfc,
        namelists,
        domain,
        zboundaries,
        +,
    )
    if ko == 0
        for k in 1:nbz
            bvsstrattfc[:, :, k] .=
                g_ndim .* p[:, :, k0 - 1] ./
                (rho[:, :, k0 - 1] .+ rhostrattfc[:, :, k0 - 1]) ./
                (thetastrattfc[:, :, k0 - 1] .^ 2) ./ jac[:, :, k0 - 1] .*
                (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
        end
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
