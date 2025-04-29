function correct!(
    state::State,
    dt::AbstractFloat,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    correct!(state, dt, U(), facray, facprs)
    correct!(state, dt, V(), facray, facprs)
    correct!(state, dt, W(), facray, facprs)
    correct!(state, dt, RhoP(), facray, facprs)
    correct!(state, PiP())
    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; model, zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; kr_sp_tfc) = state.sponge
    (; corx) = state.poisson.correction
    (; dpip) = state.variables.tendencies
    (; rho, u, p) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        facu = 1.0

        if spongelayer && sponge_uv
            facu +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k]) *
                facray
        end

        # Compute values at cell edges.
        rhou =
            0.5 * (
                rho[i, j, k] +
                rho[i + 1, j, k] +
                rhostrattfc[i, j, k] +
                rhostrattfc[i + 1, j, k]
            )
        pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
        met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

        # Compute pressure difference gradient component.
        if k == k0 && zboundaries == SolidWallBoundaries()
            dpuuedger = 0.5 * (dpip[i, j, k + 2] + dpip[i + 1, j, k + 2])
            dpuedger = 0.5 * (dpip[i, j, k + 1] + dpip[i + 1, j, k + 1])
            dpedger = 0.5 * (dpip[i, j, k] + dpip[i + 1, j, k])
            pgradx =
                kappainv * mainv2 / rhou *
                pedger *
                (
                    (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
                    met13edger *
                    (-dpuuedger + 4.0 * dpuedger - 3.0 * dpedger) *
                    0.5 / dz
                )
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            dpddedger = 0.5 * (dpip[i, j, k - 2] + dpip[i + 1, j, k - 2])
            dpdedger = 0.5 * (dpip[i, j, k - 1] + dpip[i + 1, j, k - 1])
            dpedger = 0.5 * (dpip[i, j, k] + dpip[i + 1, j, k])
            pgradx =
                kappainv * mainv2 / rhou *
                pedger *
                (
                    (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
                    met13edger *
                    (dpddedger - 4.0 * dpdedger + 3.0 * dpedger) *
                    0.5 / dz
                )
        else
            dpuedger = 0.5 * (dpip[i, j, k + 1] + dpip[i + 1, j, k + 1])
            dpdedger = 0.5 * (dpip[i, j, k - 1] + dpip[i + 1, j, k - 1])
            pgradx =
                kappainv * mainv2 / rhou *
                pedger *
                (
                    (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
                    met13edger * (dpuedger - dpdedger) * 0.5 / dz
                )
        end

        # Compute velocity correction.
        corx[i, j, k] = facprs * dt / facu * pgradx
        if model == Compressible()
            jpr =
                (
                    jac[i, j, k] * p[i, j, k] +
                    jac[i + 1, j, k] * p[i + 1, j, k]
                ) / 2
            du = -jpr * corx[i, j, k]
        else
            du = -corx[i, j, k]
        end

        u[i, j, k] += du
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; model, zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; kr_sp_tfc) = state.sponge
    (; cory) = state.poisson.correction
    (; dpip) = state.variables.tendencies
    (; rho, v, p) = state.variables.predictands

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        facv = 1.0

        if spongelayer && sponge_uv
            facv +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) *
                facray
        end

        # Compute values at cell edges.
        rhov =
            0.5 * (
                rho[i, j, k] +
                rho[i, j + 1, k] +
                rhostrattfc[i, j, k] +
                rhostrattfc[i, j + 1, k]
            )
        pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
        met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

        # Compute pressure difference gradient component.
        if k == k0 && zboundaries == SolidWallBoundaries()
            dpuuedgef = 0.5 * (dpip[i, j, k + 2] + dpip[i, j + 1, k + 2])
            dpuedgef = 0.5 * (dpip[i, j, k + 1] + dpip[i, j + 1, k + 1])
            dpedgef = 0.5 * (dpip[i, j, k] + dpip[i, j + 1, k])
            pgrady =
                kappainv * mainv2 / rhov *
                pedgef *
                (
                    (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
                    met23edgef *
                    (-dpuuedgef + 4.0 * dpuedgef - 3.0 * dpedgef) *
                    0.5 / dz
                )
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            dpddedgef = 0.5 * (dpip[i, j, k - 2] + dpip[i, j + 1, k - 2])
            dpdedgef = 0.5 * (dpip[i, j, k - 1] + dpip[i, j + 1, k - 1])
            dpedgef = 0.5 * (dpip[i, j, k] + dpip[i, j + 1, k])
            pgrady =
                kappainv * mainv2 / rhov *
                pedgef *
                (
                    (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
                    met23edgef *
                    (dpddedgef - 4.0 * dpdedgef + 3.0 * dpedgef) *
                    0.5 / dz
                )
        else
            dpuedgef = 0.5 * (dpip[i, j, k + 1] + dpip[i, j + 1, k + 1])
            dpdedgef = 0.5 * (dpip[i, j, k - 1] + dpip[i, j + 1, k - 1])
            pgrady =
                kappainv * mainv2 / rhov *
                pedgef *
                (
                    (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
                    met23edgef * (dpuedgef - dpdedgef) * 0.5 / dz
                )
        end

        # Compute velocity correction.
        cory[i, j, k] = facprs * dt / facv * pgrady
        if model == Compressible()
            jpf =
                (
                    jac[i, j, k] * p[i, j, k] +
                    jac[i, j + 1, k] * p[i, j + 1, k]
                ) / 2
            dv = -jpf * cory[i, j, k]
        else
            dv = -cory[i, j, k]
        end

        v[i, j, k] += dv
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    (; spongelayer) = state.namelists.sponge
    (; model, zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
    (; kr_sp_w_tfc) = state.sponge
    (; corx, cory) = state.poisson.correction
    (; dpip) = state.variables.tendencies
    (; rho, w, p) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in correct!: Unknown zboundaries!")
    end

    for k in k0:(k1 - 1), j in j0:j1, i in i0:i1
        facw = 1.0

        if spongelayer
            facw +=
                dt * (
                    jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] +
                    jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1]) * facray
        end

        # Compute values at cell edges.
        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhoedge =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1]) + rhostratedgeu
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        bvsstw =
            (
                jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
                jac[i, j, k] * bvsstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met13edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met33edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k + 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        dpredgeu =
            (
                jac[i + 1, j, k + 1] * dpip[i + 1, j, k] +
                jac[i + 1, j, k] * dpip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        dpledgeu =
            (
                jac[i - 1, j, k + 1] * dpip[i - 1, j, k] +
                jac[i - 1, j, k] * dpip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        dpfedgeu =
            (
                jac[i, j + 1, k + 1] * dpip[i, j + 1, k] +
                jac[i, j + 1, k] * dpip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        dpbedgeu =
            (
                jac[i, j - 1, k + 1] * dpip[i, j - 1, k] +
                jac[i, j - 1, k] * dpip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

        # Compute pressure difference gradient component.
        pgradz =
            kappainv * mainv2 / rhoedge *
            pedgeu *
            (
                met13edgeu * (dpredgeu - dpledgeu) * 0.5 / dx +
                met23edgeu * (dpfedgeu - dpbedgeu) * 0.5 / dy +
                met33edgeu * (dpip[i, j, k + 1] - dpip[i, j, k]) / dz
            )

        # Compute velocity correction.
        if model == Compressible()
            jpu =
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (p[i, j, k] + p[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])

            dw =
                -facprs * dt / (facw + bvsstw * dt^2.0) * jpu * pgradz -
                1.0 / (facw + bvsstw * dt^2.0) *
                bvsstw *
                dt^2.0 *
                jpu *
                0.5 *
                (
                    jac[i, j, k + 1] * (
                        met[i, j, k, 1, 3] *
                        (corx[i, j, k] + corx[i - 1, j, k]) +
                        met[i, j, k, 2, 3] *
                        (cory[i, j, k] + cory[i, j - 1, k])
                    ) +
                    jac[i, j, k] * (
                        met[i, j, k + 1, 1, 3] *
                        (corx[i, j, k + 1] + corx[i - 1, j, k + 1]) +
                        met[i, j, k + 1, 2, 3] *
                        (cory[i, j, k + 1] + cory[i, j - 1, k + 1])
                    )
                ) / (jac[i, j, k] + jac[i, j, k + 1])
        else
            dw =
                -facprs * dt /
                (facw + rhostratedgeu / rhoedge * bvsstw * dt^2.0) * pgradz -
                1.0 / (facw + rhostratedgeu / rhoedge * bvsstw * dt^2.0) *
                rhostratedgeu / rhoedge *
                bvsstw *
                dt^2.0 *
                0.5 *
                (
                    jac[i, j, k + 1] * (
                        met[i, j, k, 1, 3] *
                        (corx[i, j, k] + corx[i - 1, j, k]) +
                        met[i, j, k, 2, 3] *
                        (cory[i, j, k] + cory[i, j - 1, k])
                    ) +
                    jac[i, j, k] * (
                        met[i, j, k + 1, 1, 3] *
                        (corx[i, j, k + 1] + corx[i - 1, j, k + 1]) +
                        met[i, j, k + 1, 2, 3] *
                        (cory[i, j, k + 1] + cory[i, j - 1, k + 1])
                    )
                ) / (jac[i, j, k] + jac[i, j, k + 1])
        end

        w[i, j, k] += dw
    end

    return
end

function correct!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    (; spongelayer) = state.namelists.sponge
    (; model, zboundaries) = state.namelists.setting
    (; kappainv, mainv2, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
    (; kr_sp_w_tfc) = state.sponge
    (; corx, cory) = state.poisson.correction
    (; dpip) = state.variables.tendencies
    (; rho, rhop) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        facw = 1.0

        if spongelayer
            facw += dt * kr_sp_w_tfc[i, j, k] * facray
        end

        # Compute P coefficients.
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        pedged =
            (
                jac[i, j, k - 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Compute density coefficients.
        rhow0 =
            (
                jac[i, j, k + 1] * (rho[i, j, k] + rhostrattfc[i, j, k]) +
                jac[i, j, k] * (rho[i, j, k + 1] + rhostrattfc[i, j, k + 1])
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhowm =
            (
                jac[i, j, k - 1] * (rho[i, j, k] + rhostrattfc[i, j, k]) +
                jac[i, j, k] * (rho[i, j, k - 1] + rhostrattfc[i, j, k - 1])
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Interpolate metric tensor elements.
        met13edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met13edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k - 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met23edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k - 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met33edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k + 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met33edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k - 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Interpolate pressure differences.
        dpredgeu =
            (
                jac[i + 1, j, k + 1] * dpip[i + 1, j, k] +
                jac[i + 1, j, k] * dpip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        dpledgeu =
            (
                jac[i - 1, j, k + 1] * dpip[i - 1, j, k] +
                jac[i - 1, j, k] * dpip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        dpredged =
            (
                jac[i + 1, j, k - 1] * dpip[i + 1, j, k] +
                jac[i + 1, j, k] * dpip[i + 1, j, k - 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
        dpledged =
            (
                jac[i - 1, j, k - 1] * dpip[i - 1, j, k] +
                jac[i - 1, j, k] * dpip[i - 1, j, k - 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
        dpfedgeu =
            (
                jac[i, j + 1, k + 1] * dpip[i, j + 1, k] +
                jac[i, j + 1, k] * dpip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        dpbedgeu =
            (
                jac[i, j - 1, k + 1] * dpip[i, j - 1, k] +
                jac[i, j - 1, k] * dpip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
        dpfedged =
            (
                jac[i, j + 1, k - 1] * dpip[i, j + 1, k] +
                jac[i, j + 1, k] * dpip[i, j + 1, k - 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
        dpbedged =
            (
                jac[i, j - 1, k - 1] * dpip[i, j - 1, k] +
                jac[i, j - 1, k] * dpip[i, j - 1, k - 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])

        # Compute pressure difference gradients.
        pgradzedgeu =
            kappainv * mainv2 * pedgeu / rhow0 * (
                0.5 * met13edgeu * (dpredgeu - dpledgeu) / dx +
                0.5 * met23edgeu * (dpfedgeu - dpbedgeu) / dy +
                met33edgeu * (dpip[i, j, k + 1] - dpip[i, j, k]) / dz
            )
        pgradzedged =
            kappainv * mainv2 * pedged / rhowm * (
                0.5 * met13edged * (dpredged - dpledged) / dx +
                0.5 * met23edged * (dpfedged - dpbedged) / dy +
                met33edged * (dpip[i, j, k] - dpip[i, j, k - 1]) / dz
            )

        # Adjust at boundaries.
        if k == k0 && zboundaries == SolidWallBoundaries()
            pgradzedged = 0.0
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            pgradzedgeu = 0.0
        end

        # Interpolate.
        pgradz = 0.5 * (pgradzedgeu + pgradzedged)

        # Compute buoyancy correction.
        if model == Compressible()
            db =
                -1.0 / (facw + bvsstrattfc[i, j, k] * dt^2.0) * (
                    -bvsstrattfc[i, j, k] *
                    facprs *
                    dt^2.0 *
                    jac[i, j, k] *
                    pgradz +
                    bvsstrattfc[i, j, k] *
                    dt *
                    jac[i, j, k] *
                    facw *
                    0.5 *
                    (
                        met[i, j, k, 1, 3] *
                        (corx[i, j, k] + corx[i - 1, j, k]) +
                        met[i, j, k, 2, 3] *
                        (cory[i, j, k] + cory[i, j - 1, k])
                    )
                )
        else
            db =
                -1.0 / (
                    facw +
                    rhostrattfc[i, j, k] /
                    (rho[i, j, k] + rhostrattfc[i, j, k]) *
                    bvsstrattfc[i, j, k] *
                    dt^2.0
                ) * (
                    -rhostrattfc[i, j, k] /
                    (rho[i, j, k] + rhostrattfc[i, j, k]) *
                    bvsstrattfc[i, j, k] *
                    facprs *
                    dt^2.0 *
                    jac[i, j, k] *
                    pgradz +
                    rhostrattfc[i, j, k] /
                    (rho[i, j, k] + rhostrattfc[i, j, k]) *
                    bvsstrattfc[i, j, k] *
                    dt *
                    jac[i, j, k] *
                    facw *
                    0.5 *
                    (
                        met[i, j, k, 1, 3] *
                        (corx[i, j, k] + corx[i - 1, j, k]) +
                        met[i, j, k, 2, 3] *
                        (cory[i, j, k] + cory[i, j - 1, k])
                    )
                )
        end

        rhop[i, j, k] -= (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim * db
    end

    return
end

function correct!(state::State, variable::PiP)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; pip) = state.variables.predictands
    (; dpip) = state.variables.tendencies

    @views pip[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .=
        pip[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .+
        dpip[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]

    return
end
