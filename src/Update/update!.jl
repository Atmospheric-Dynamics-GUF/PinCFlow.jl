function update!(state::State, dt::AbstractFloat, m::Integer, variable::Rho)
    (; model) = state.namelists.setting
    update!(state, dt, m, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Boussinesq,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::PseudoIncompressible,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; drho) = state.variables.tendencies
    (; phirho) = state.variables.fluxes
    (; rho) = state.variables.predictands

    if m == 1
        drho .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phirho[i - 1, j, k, 1]
        fr = phirho[i, j, k, 1]
        gb = phirho[i, j - 1, k, 2]
        gf = phirho[i, j, k, 2]
        hd = phirho[i, j, k - 1, 3]
        hu = phirho[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        f = -fluxdiff

        drho[i, j, k] = dt * f + alphark[m] * drho[i, j, k]
        rho[i, j, k] += betark[m] * drho[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::RhoP,
    side::LHS,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; drhop) = state.variables.tendencies
    (; phirhop) = state.variables.fluxes
    (; rhop) = state.variables.predictands

    if m == 1
        drhop .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phirhop[i - 1, j, k, 1]
        fr = phirhop[i, j, k, 1]
        gb = phirhop[i, j - 1, k, 2]
        gf = phirhop[i, j, k, 2]
        hd = phirhop[i, j, k - 1, 3]
        hu = phirhop[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        f = -fluxdiff

        drhop[i, j, k] = dt * f + alphark[m] * drhop[i, j, k]
        rhop[i, j, k] += betark[m] * drhop[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::EXPL,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; grid) = state
    (; g_ndim) = state.constants
    (; rhostrattfc, bvsstrattfc) = state.atmosphere
    (; predictands) = state.variables
    (; rho, rhop) = predictands

    for k in k0:k1, j in j0:j1, i in i0:i1
        wvrt =
            0.5 * (
                compute_vertical_wind(i, j, k, predictands, grid) +
                compute_vertical_wind(i, j, k - 1, predictands, grid)
            )
        buoy = -g_ndim * rhop[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k])
        buoy -=
            dt * rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
            bvsstrattfc[i, j, k] *
            wvrt
        rhop[i, j, k] = -buoy * (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::IMPL,
    facray::AbstractFloat,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; zboundaries) = state.namelists.setting
    (; spongelayer) = state.namelists.sponge
    (; kr_sp_w_tfc) = state.sponge
    (; kappainv, mainv2, g_ndim) = state.constants
    (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
    (; rho, rhop, u, v, pip) = state.variables.predictands
    (; wold) = state.variables.backups

    for k in k0:k1, j in j0:j1, i in i0:i1
        rhow =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhowm =
            (
                jac[i, j, k - 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        rhow +=
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhowm +=
            (
                jac[i, j, k - 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Momentum is predicted before buoyancy in implicit
        # steps.
        wvrt = 0.5 * (wold[i, j, k] + wold[i, j, k - 1])

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

        # Interpolate metric-tensor elements.
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
        met13edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k - 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met23edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k - 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met33edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k - 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        # Interpolate pressure differences.
        piredgeu =
            (
                jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        piledgeu =
            (
                jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        piredged =
            (
                jac[i + 1, j, k - 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k - 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
        piledged =
            (
                jac[i - 1, j, k - 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k - 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
        pifedgeu =
            (
                jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        pibedgeu =
            (
                jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
        pifedged =
            (
                jac[i, j + 1, k - 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k - 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
        pibedged =
            (
                jac[i, j - 1, k - 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k - 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])

        # Compute pressure gradients.
        pigradzedgeu =
            kappainv * mainv2 * pedgeu / rhow * (
                0.5 * met13edgeu * (piredgeu - piledgeu) / dx +
                0.5 * met23edgeu * (pifedgeu - pibedgeu) / dy +
                met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
            )
        pigradzedged =
            kappainv * mainv2 * pedged / rhowm * (
                0.5 * met13edged * (piredged - piledged) / dx +
                0.5 * met23edged * (pifedged - pibedged) / dy +
                met33edged * (pip[i, j, k] - pip[i, j, k - 1]) / dz
            )

        # Adjust at boundaries.
        if k == k0 && zboundaries == SolidWallBoundaries()
            pigradzedged = 0.0
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            pigradzedgeu = 0.0
        end

        # Interpolate.
        pigrad = 0.5 * (pigradzedgeu + pigradzedged)

        facw = 1.0

        if spongelayer
            facw += dt * kr_sp_w_tfc[i, j, k] * facray
        end

        # Predict buoyancy.
        buoy = -g_ndim * rhop[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k])
        buoy =
            1.0 / (
                facw +
                rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
                bvsstrattfc[i, j, k] *
                dt^2.0
            ) * (
                -rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
                bvsstrattfc[i, j, k] *
                dt *
                jac[i, j, k] *
                (wvrt - dt * pigrad) +
                facw * buoy +
                rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
                bvsstrattfc[i, j, k] *
                dt *
                jac[i, j, k] *
                facw *
                0.5 *
                (
                    met[i, j, k, 1, 3] * (u[i, j, k] + u[i - 1, j, k]) +
                    met[i, j, k, 2, 3] * (v[i, j, k] + v[i, j - 1, k])
                )
            )

        rhop[i, j, k] = -buoy * (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::U,
    side::LHS,
)
    (; alphark, betark) = state.time
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhostrattfc, fc) = state.atmosphere
    (; du) = state.variables.tendencies
    (; phiu) = state.variables.fluxes
    (; rhoold, uold) = state.variables.backups
    (; rho, u, v) = state.variables.predictands

    if m == 1
        du .= 0.0
    end

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1

        # Compute zonal momentum flux divergence.
        fr = phiu[i, j, k, 1]
        fl = phiu[i - 1, j, k, 1]
        gf = phiu[i, j, k, 2]
        gb = phiu[i, j - 1, k, 2]
        hu = phiu[i, j, k, 3]
        hd = phiu[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust zonal momentum flux divergence.
        jacedger = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])
        fluxdiff /= jacedger

        # Explicit integration of Coriolis force in TFC.
        uold[i, j, k] = u[i, j, k]
        vc = 0.5 * (v[i, j, k] + v[i, j - 1, k])
        vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
        volforce =
            0.5 *
            fc[j] *
            (
                (rhoold[i, j, k] + rhostrattfc[i, j, k]) * vc +
                (rhoold[i + 1, j, k] + rhostrattfc[i + 1, j, k]) * vr
            )

        # Compute force.
        force = -fluxdiff + volforce

        # Interpolate density.
        rhom_1 = 0.5 * (rhoold[i, j, k] + rhoold[i + 1, j, k])
        rhom = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        rhom_1 += rhostratedger
        rhom += rhostratedger

        # Set velocity and momentum at previous time.
        um_1 = u[i, j, k]
        momm_1 = rhom_1 * um_1

        # Compute tendency.
        du[i, j, k] = dt * force + alphark[m] * du[i, j, k]

        # Update momentum.
        momm = momm_1 + betark[m] * du[i, j, k]

        # Update wind.
        uast = momm / rhom
        u[i, j, k] = uast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::EXPL,
)
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho, u, pip) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhou = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        rhou += rhostratedger

        pir = pip[i + 1, j, k]
        pil = pip[i, j, k]

        # Compute values at cell edges.
        pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
        met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

        # Compute pressure gradient component.
        if k == k0 && zboundaries == SolidWallBoundaries()
            piuuedger = 0.5 * (pip[i, j, k + 2] + pip[i + 1, j, k + 2])
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigrad =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (-piuuedger + 4.0 * piuedger - 3.0 * piedger) *
                    0.5 / dz
                )
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            piddedger = 0.5 * (pip[i, j, k - 2] + pip[i + 1, j, k - 2])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigrad =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (piddedger - 4.0 * pidedger + 3.0 * piedger) *
                    0.5 / dz
                )
        else
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            pigrad =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger * (piuedger - pidedger) * 0.5 / dz
                )
        end

        volfcx = compute_volume_force(state, (i, j, k), U())

        uhorx = u[i, j, k]

        # Update wind.
        uast = uhorx + dt * (-pigrad + volfcx / rhou)
        u[i, j, k] = uast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::IMPL,
    facray::AbstractFloat,
)
    (; zboundaries) = state.namelists.setting
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; kr_sp_tfc) = state.sponge
    (; rho, u, pip) = state.variables.predictands

    for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhou = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhostratedger = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i + 1, j, k])
        rhou += rhostratedger

        pir = pip[i + 1, j, k]
        pil = pip[i, j, k]

        # Compute values at cell edges.
        pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
        met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

        # Compute pressure gradient component.
        if k == k0 && zboundaries == SolidWallBoundaries()
            piuuedger = 0.5 * (pip[i, j, k + 2] + pip[i + 1, j, k + 2])
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigradx =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (-piuuedger + 4.0 * piuedger - 3.0 * piedger) *
                    0.5 / dz
                )
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            piddedger = 0.5 * (pip[i, j, k - 2] + pip[i + 1, j, k - 2])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            piedger = 0.5 * (pip[i, j, k] + pip[i + 1, j, k])
            pigradx =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger *
                    (piddedger - 4.0 * pidedger + 3.0 * piedger) *
                    0.5 / dz
                )
        else
            piuedger = 0.5 * (pip[i, j, k + 1] + pip[i + 1, j, k + 1])
            pidedger = 0.5 * (pip[i, j, k - 1] + pip[i + 1, j, k - 1])
            pigradx =
                kappainv * mainv2 * pedger / rhou * (
                    (pir - pil) / dx +
                    met13edger * (piuedger - pidedger) * 0.5 / dz
                )
        end

        volfcx = compute_volume_force(state, (i, j, k), U())

        uhorx = u[i, j, k]

        facu = 1.0

        if spongelayer && sponge_uv
            facu +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k]) *
                facray
        end

        # Update wind.
        uast = 1.0 / facu * (uhorx + dt * (-pigradx + volfcx / rhou))
        u[i, j, k] = uast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::V,
    side::LHS,
)
    (; alphark, betark) = state.time
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhostrattfc, fc) = state.atmosphere
    (; dv) = state.variables.tendencies
    (; phiv) = state.variables.fluxes
    (; rhoold, uold, vold) = state.variables.backups
    (; rho, v) = state.variables.predictands

    if m == 1
        dv .= 0.0
    end

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1

        # Compute meridional momentum flux divergence.
        fr = phiv[i, j, k, 1]
        fl = phiv[i - 1, j, k, 1]
        gf = phiv[i, j, k, 2]
        gb = phiv[i, j - 1, k, 2]
        hu = phiv[i, j, k, 3]
        hd = phiv[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust meridional momentum flux divergence.
        jacedgef = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])
        fluxdiff /= jacedgef

        # Explicit integration of Coriolis force in TFC.
        vold[i, j, k] = v[i, j, k]
        uc = 0.5 * (uold[i, j, k] + uold[i - 1, j, k])
        uf = 0.5 * (uold[i, j + 1, k] + uold[i - 1, j + 1, k])

        volforce =
            -0.5 * (
                fc[j] * (rhoold[i, j, k] + rhostrattfc[i, j, k]) * uc +
                fc[j + 1] *
                (rhoold[i, j + 1, k] + rhostrattfc[i, j + 1, k]) *
                uf
            )

        force = -fluxdiff + volforce

        # Interpolate density.
        rhom_1 = 0.5 * (rhoold[i, j, k] + rhoold[i, j + 1, k])
        rhom = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        rhom_1 += rhostratedgef
        rhom += rhostratedgef

        vm_1 = v[i, j, k]
        momm_1 = rhom_1 * vm_1

        dv[i, j, k] = dt * force + alphark[m] * dv[i, j, k]

        momm = momm_1 + betark[m] * dv[i, j, k]

        # Update wind.
        vast = momm / rhom
        v[i, j, k] = vast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::EXPL,
)
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dy, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rho, v, pip) = state.variables.predictands

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhov = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        rhov += rhostratedgef

        pif = pip[i, j + 1, k]
        pib = pip[i, j, k]

        # Compute values at cell edges.
        pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
        met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

        # Compute pressure gradient component.
        if k == k0 && zboundaries == SolidWallBoundaries()
            piuuedgef = 0.5 * (pip[i, j, k + 2] + pip[i, j + 1, k + 2])
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrad =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (-piuuedgef + 4.0 * piuedgef - 3.0 * piedgef) *
                    0.5 / dz
                )
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            piddedgef = 0.5 * (pip[i, j, k - 2] + pip[i, j + 1, k - 2])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrad =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (piddedgef - 4.0 * pidedgef + 3.0 * piedgef) *
                    0.5 / dz
                )
        else
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            pigrad =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef * (piuedgef - pidedgef) * 0.5 / dz
                )
        end

        volfcy = compute_volume_force(state, (i, j, k), V())

        vhory = v[i, j, k]

        # Update wind.
        vast = vhory + dt * (-pigrad + volfcy / rhov)
        v[i, j, k] = vast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::IMPL,
    facray::AbstractFloat,
)
    (; zboundaries) = state.namelists.setting
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; kappainv, mainv2) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dy, dz, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; kr_sp_tfc) = state.sponge
    (; rho, v, pip) = state.variables.predictands

    for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhov = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhostratedgef = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j + 1, k])
        rhov += rhostratedgef

        pif = pip[i, j + 1, k]
        pib = pip[i, j, k]

        # Compute values at cell edges.
        pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
        met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

        # Compute pressure gradient component.
        if k == k0 && zboundaries == SolidWallBoundaries()
            piuuedgef = 0.5 * (pip[i, j, k + 2] + pip[i, j + 1, k + 2])
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrady =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (-piuuedgef + 4.0 * piuedgef - 3.0 * piedgef) *
                    0.5 / dz
                )
        elseif k == k1 && zboundaries == SolidWallBoundaries()
            piddedgef = 0.5 * (pip[i, j, k - 2] + pip[i, j + 1, k - 2])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            piedgef = 0.5 * (pip[i, j, k] + pip[i, j + 1, k])
            pigrady =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef *
                    (piddedgef - 4.0 * pidedgef + 3.0 * piedgef) *
                    0.5 / dz
                )
        else
            piuedgef = 0.5 * (pip[i, j, k + 1] + pip[i, j + 1, k + 1])
            pidedgef = 0.5 * (pip[i, j, k - 1] + pip[i, j + 1, k - 1])
            pigrady =
                kappainv * mainv2 * pedgef / rhov * (
                    (pif - pib) / dy +
                    met23edgef * (piuedgef - pidedgef) * 0.5 / dz
                )
        end

        volfcy = compute_volume_force(state, (i, j, k), V())

        vhory = v[i, j, k]

        facv = 1.0

        if spongelayer && sponge_uv
            facv +=
                dt *
                0.5 *
                (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) *
                facray
        end

        # Update wind.
        vast = 1.0 / facv * (vhory + dt * (-pigrady + volfcy / rhov))
        v[i, j, k] = vast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::W,
    side::LHS,
)
    (; zboundaries) = state.namelists.setting
    (; alphark, betark) = state.time
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; grid) = state
    (; dx, dy, dz, jac, met) = grid
    (; rhostrattfc, fc) = state.atmosphere
    (; dw) = state.variables.tendencies
    (; phiu, phiv, phiw) = state.variables.fluxes
    (; rhoold, uold, vold) = state.variables.backups
    (; rho, w) = state.variables.predictands

    # Initialize fields for transformation of momentum flux divergence.
    (fluxdiffu, fluxdiffv) = (zeros(2, 2) for i in 1:2)

    if m == 1
        dw .= 0.0
    end

    if zboundaries != SolidWallBoundaries()
        error("Error in update_momentum!: Unknown case zBoundary!")
    end

    for k in k0:(k1 - 1), j in j0:j1, i in i0:i1
        # Compute vertical momentum flux divergence.
        fr = phiw[i, j, k, 1]
        fl = phiw[i - 1, j, k, 1]
        gf = phiw[i, j, k, 2]
        gb = phiw[i, j - 1, k, 2]
        hu = phiw[i, j, k, 3]
        hd = phiw[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust Cartesian vertical momentum flux divergence.
        jacedgeu =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])
        fluxdiff /= jacedgeu

        # Compute zonal momentum flux divergences.
        for ll in 0:1, mm in 0:1
            fr = phiu[i - ll, j, k + mm, 1]
            fl = phiu[i - 1 - ll, j, k + mm, 1]
            gf = phiu[i - ll, j, k + mm, 2]
            gb = phiu[i - ll, j - 1, k + mm, 2]
            hu = phiu[i - ll, j, k + mm, 3]
            hd = phiu[i - ll, j, k - 1 + mm, 3]
            fluxdiffu[ll + 1, mm + 1] =
                (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            jacedger =
                0.5 * (jac[i - ll, j, k + mm] + jac[i + 1 - ll, j, k + mm])
            fluxdiffu[ll + 1, mm + 1] /= jacedger
        end

        # Compute meridional momentum flux divergences.
        for ll in 0:1, mm in 0:1
            fr = phiv[i, j - ll, k + mm, 1]
            fl = phiv[i - 1, j - ll, k + mm, 1]
            gf = phiv[i, j - ll, k + mm, 2]
            gb = phiv[i, j - 1 - ll, k + mm, 2]
            hu = phiv[i, j - ll, k + mm, 3]
            hd = phiv[i, j - ll, k - 1 + mm, 3]
            fluxdiffv[ll + 1, mm + 1] =
                (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            jacedgef =
                0.5 * (jac[i, j - ll, k + mm] + jac[i, j + 1 - ll, k + mm])
            fluxdiffv[ll + 1, mm + 1] /= jacedgef
        end

        # Compute transformed vertical momentum flux divergence.
        fluxdiff = transform(
            i,
            j,
            k,
            fluxdiffu[1, 1],
            fluxdiffu[1, 2],
            fluxdiffu[2, 1],
            fluxdiffu[2, 2],
            fluxdiffv[1, 1],
            fluxdiffv[1, 2],
            fluxdiffv[2, 1],
            fluxdiffv[2, 2],
            fluxdiff,
            TFC(),
            grid,
        )

        # Explicit integration of Coriolis force in TFC.
        vc = 0.5 * (vold[i, j, k] + vold[i, j - 1, k])
        vu = 0.5 * (vold[i, j, k + 1] + vold[i, j - 1, k + 1])
        uc = 0.5 * (uold[i, j, k] + uold[i - 1, j, k])
        uu = 0.5 * (uold[i, j, k + 1] + uold[i - 1, j, k + 1])

        volforce =
            fc[j] * (
                jac[i, j, k + 1] *
                met[i, j, k, 1, 3] *
                (rhoold[i, j, k] + rhostrattfc[i, j, k]) *
                vc +
                jac[i, j, k] *
                met[i, j, k + 1, 1, 3] *
                (rhoold[i, j, k + 1] + rhostrattfc[i, j, k + 1]) *
                vu
            ) / (jac[i, j, k] + jac[i, j, k + 1]) -
            fc[j] * (
                jac[i, j, k + 1] *
                met[i, j, k, 2, 3] *
                (rhoold[i, j, k] + rhostrattfc[i, j, k]) *
                uc +
                jac[i, j, k] *
                met[i, j, k + 1, 2, 3] *
                (rhoold[i, j, k + 1] + rhostrattfc[i, j, k + 1]) *
                uu
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        force = -fluxdiff + volforce

        # Interpolate densities.
        rhom_1 =
            (
                jac[i, j, k + 1] * rhoold[i, j, k] +
                jac[i, j, k] * rhoold[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhom =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhom_1 += rhostratedgeu
        rhom += rhostratedgeu

        wm_1 = w[i, j, k]
        momm_1 = rhom_1 * wm_1

        dw[i, j, k] = dt * force + alphark[m] * dw[i, j, k]

        momm = momm_1 + betark[m] * dw[i, j, k]

        # Update wind.
        wast = momm / rhom
        w[i, j, k] = wast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::EXPL,
)
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; rhopold) = state.variables.backups
    (; rho, w, pip) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in update!: Unknown zboundaries!")
    end

    for k in k0:(k1 - 1), j in j0:j1, i in i0:i1
        rho000 = rho[i, j, k]
        rho001 = rho[i, j, k + 1]

        rhow =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rho000 += rhostrattfc[i, j, k]
        rho001 += rhostrattfc[i, j, k + 1]
        rhow += rhostratedgeu

        # Compute values at cell edges.
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
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
        piredgeu =
            (
                jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        piledgeu =
            (
                jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        pifedgeu =
            (
                jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        pibedgeu =
            (
                jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

        # Compute pressure gradient component.
        pigrad =
            kappainv * mainv2 * pedgeu / rhow * (
                met13edgeu * (piredgeu - piledgeu) * 0.5 / dx +
                met23edgeu * (pifedgeu - pibedgeu) * 0.5 / dy +
                met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
            )

        volfcz = compute_volume_force(state, (i, j, k), W())

        wvert = w[i, j, k]

        buoy =
            -g_ndim * (
                jac[i, j, k + 1] * rhopold[i, j, k] / rho000 / jac[i, j, k] +
                jac[i, j, k] * rhopold[i, j, k + 1] / rho001 / jac[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        # Update wind.
        wast = wvert + dt * (buoy - pigrad + volfcz / rhow)
        w[i, j, k] = wast
    end

    # Return.
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::IMPL,
    facray::AbstractFloat,
)
    (; spongelayer) = state.namelists.sponge
    (; zboundaries) = state.namelists.setting
    (; kappainv, mainv2, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
    (; kr_sp_w_tfc) = state.sponge
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in update!: Unknown zboundaries!")
    end

    for k in k0:(k1 - 1), j in j0:j1, i in i0:i1
        rho000 = rho[i, j, k]
        rho001 = rho[i, j, k + 1]

        rhow =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        rhostratedgeu =
            (
                jac[i, j, k + 1] * rhostrattfc[i, j, k] +
                jac[i, j, k] * rhostrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rho000 += rhostrattfc[i, j, k]
        rho001 += rhostrattfc[i, j, k + 1]
        rhow += rhostratedgeu

        # Compute values at cell edges.
        pedgeu =
            (
                jac[i, j, k + 1] * pstrattfc[i, j, k] +
                jac[i, j, k] * pstrattfc[i, j, k + 1]
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
        piredgeu =
            (
                jac[i + 1, j, k + 1] * pip[i + 1, j, k] +
                jac[i + 1, j, k] * pip[i + 1, j, k + 1]
            ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        piledgeu =
            (
                jac[i - 1, j, k + 1] * pip[i - 1, j, k] +
                jac[i - 1, j, k] * pip[i - 1, j, k + 1]
            ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
        pifedgeu =
            (
                jac[i, j + 1, k + 1] * pip[i, j + 1, k] +
                jac[i, j + 1, k] * pip[i, j + 1, k + 1]
            ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        pibedgeu =
            (
                jac[i, j - 1, k + 1] * pip[i, j - 1, k] +
                jac[i, j - 1, k] * pip[i, j - 1, k + 1]
            ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

        # Compute pressure gradient component.
        pigrad =
            kappainv * mainv2 * pedgeu / rhow * (
                met13edgeu * (piredgeu - piledgeu) * 0.5 / dx +
                met23edgeu * (pifedgeu - pibedgeu) * 0.5 / dy +
                met33edgeu * (pip[i, j, k + 1] - pip[i, j, k]) / dz
            )

        volfcz = compute_volume_force(state, (i, j, k), W())

        wvert = w[i, j, k]

        bvsstw =
            (
                jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
                jac[i, j, k] * bvsstrattfc[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        facw = 1.0

        if spongelayer
            facw +=
                dt * (
                    jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] +
                    jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1]) * facray
        end

        # Buoyancy is predicted after momentum in implicit steps.
        buoy =
            -g_ndim * (
                jac[i, j, k + 1] * rhop[i, j, k] / rho000 / jac[i, j, k] +
                jac[i, j, k] * rhop[i, j, k + 1] / rho001 / jac[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        uc = 0.5 * (u[i, j, k] + u[i - 1, j, k])
        uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
        vc = 0.5 * (v[i, j, k] + v[i, j - 1, k])
        vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])

        # Update wind.
        wast =
            1.0 / (facw + rhostratedgeu / rhow * bvsstw * dt^2.0) * (
                wvert - dt * pigrad +
                dt * buoy +
                dt * volfcz / rhow +
                rhostratedgeu / rhow *
                bvsstw *
                dt^2.0 *
                (
                    jac[i, j, k + 1] *
                    (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
                    jac[i, j, k] *
                    (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
                ) / (jac[i, j, k] + jac[i, j, k + 1])
            )
        w[i, j, k] = wast
    end

    # Return.
    return
end
