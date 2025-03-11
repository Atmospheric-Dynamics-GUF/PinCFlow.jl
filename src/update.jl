function massUpdate!(model, ode, dt, upd_var, upd_mod, int_mod, m, facray)
    @trixi_timeit timer() "Mass update" begin
    #! format: noindent
    (; kr_sp_tfc, kr_sp_w_tfc) = model.atmosphere.sponge
    var = model.variables.prognostic_fields
    flux = model.fluxes
    (;
        pstrattfc,
        rhostrattfc,
        bvsstrattfc
    ) = model.atmosphere
    jac = model.grid.jac
    met = model.grid.met

    (; nx, ny, nz) = model.domain
    (; dx, dy, dz) = model.grid
    dRho = model.variables.tendencies.drho
    dRhop = model.variables.tendencies.drhop
    wOld = model.variables.history.w
    (; alphark, betark) = ode
    (; g_ndim, kappainv, mainv2) = model.constants

    if m == 1 && upd_var == "rho"
        dRho .= 0.0
    elseif m == 1 && upd_var == "rhop"
        dRhop .= 0.0
    end

    if upd_var == "rho"
        if !(upd_mod in ["tot", "lhs"])
            error("ERROR: wrong upd_mod for upd_var = rho")
        end

        if int_mod != "expl"
            error("ERROR: wrong int_mod for upd_var = rho")
        end

        for k in 1:nz, j in 1:ny, i in 1:nx
            fL = flux.rho[i-1, j, k, 1]
            fR = flux.rho[i, j, k, 1]
            gB = flux.rho[i, j-1, k, 2]
            gF = flux.rho[i, j, k, 2]
            hD = flux.rho[i, j, k-1, 3]
            hU = flux.rho[i, j, k, 3]
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
            fluxDiff /= jac[i, j, k]
            F = -fluxDiff

            dRho[i, j, k] = dt * F + alphark[m] * dRho[i, j, k]
            var.rho[i, j, k] += betark[m] * dRho[i, j, k]
        end
    elseif upd_var == "rhop"
        if upd_mod == "lhs"
            if int_mod != "expl"
                error("ERROR: wrong int_mod for upd_mod = lhs")
            end

            for k in 1:nz, j in 1:ny, i in 1:nx
                fL = flux.rhop[i-1, j, k, 1]
                fR = flux.rhop[i, j, k, 1]
                gB = flux.rhop[i, j-1, k, 2]
                gF = flux.rhop[i, j, k, 2]
                hD = flux.rhop[i, j, k-1, 3]
                hU = flux.rhop[i, j, k, 3]

                fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
                fluxDiff /= jac[i, j, k]
                F = -fluxDiff

                dRhop[i, j, k] = dt * F + alphark[m] * dRhop[i, j, k]
                var.rhop[i, j, k] += betark[m] * dRhop[i, j, k]
            end
        elseif upd_mod == "rhs"
            if int_mod == "impl"
                kr_sp_tfc .*= facray
                kr_sp_w_tfc .*= facray

                for k in 1:nz, j in 1:ny, i in 1:nx
                    rho = var.rho[i, j, k]
                    rhow = (jac[i, j, k+1] * var.rho[i, j, k] +
                            jac[i, j, k] * var.rho[i, j, k+1]) /
                           (jac[i, j, k] + jac[i, j, k+1])
                    rhowm = (jac[i, j, k-1] * var.rho[i, j, k] +
                             jac[i, j, k] * var.rho[i, j, k-1]) /
                            (jac[i, j, k] + jac[i, j, k-1])

                    rho += rhostrattfc[i, j, k]
                    rhow += (jac[i, j, k+1] * rhostrattfc[i, j, k] +
                             jac[i, j, k] * rhostrattfc[i, j, k+1]) /
                            (jac[i, j, k] + jac[i, j, k+1])
                    rhowm += (jac[i, j, k-1] * rhostrattfc[i, j, k] +
                              jac[i, j, k] * rhostrattfc[i, j, k-1]) /
                             (jac[i, j, k] + jac[i, j, k-1])

                    wvrt = 0.5 * (wOld[i, j, k] + wOld[i, j, k-1])

                    pEdgeU = (jac[i, j, k+1] * pstrattfc[i, j, k] +
                              jac[i, j, k] * pstrattfc[i, j, k+1]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                    pEdgeD = (jac[i, j, k-1] * pstrattfc[i, j, k] +
                              jac[i, j, k] * pstrattfc[i, j, k-1]) /
                             (jac[i, j, k] + jac[i, j, k-1])

                    # Interpolate metric-tensor elements.
                    met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                                  jac(i, j, k) * met(i, j, k + 1, 1, 3)) /
                                 (jac(i, j, k) + jac(i, j, k + 1))
                    met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                                  jac(i, j, k) * met(i, j, k + 1, 2, 3)) /
                                 (jac(i, j, k) + jac(i, j, k + 1))
                    met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                                  jac(i, j, k) * met(i, j, k + 1, 3, 3)) /
                                 (jac(i, j, k) + jac(i, j, k + 1))
                    met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) +
                                  jac(i, j, k) * met(i, j, k - 1, 1, 3)) /
                                 (jac(i, j, k) + jac(i, j, k - 1))
                    met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) +
                                  jac(i, j, k) * met(i, j, k - 1, 2, 3)) /
                                 (jac(i, j, k) + jac(i, j, k - 1))
                    met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) +
                                  jac(i, j, k) * met(i, j, k - 1, 3, 3)) /
                                 (jac(i, j, k) + jac(i, j, k - 1))

                    # Interpolate pressure differences.
                    piREdgeU = (jac(i + 1, j, k + 1) * var.pip(i + 1, j, k) +
                                jac(i + 1, j, k) * var.pip(i + 1, j, k + 1)) /
                               (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
                    piLEdgeU = (jac(i - 1, j, k + 1) * var.pip(i - 1, j, k) +
                                jac(i - 1, j, k) * var.pip(i - 1, j, k + 1)) /
                               (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    piREdgeD = (jac(i + 1, j, k - 1) * var.pip(i + 1, j, k) +
                                jac(i + 1, j, k) * var.pip(i + 1, j, k - 1)) /
                               (jac(i + 1, j, k) + jac(i + 1, j, k - 1))
                    piLEdgeD = (jac(i - 1, j, k - 1) * var.pip(i - 1, j, k) +
                                jac(i - 1, j, k) * var.pip(i - 1, j, k - 1)) /
                               (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
                    piFEdgeU = (jac(i, j + 1, k + 1) * var.pip(i, j + 1, k) +
                                jac(i, j + 1, k) * var.pip(i, j + 1, k + 1)) /
                               (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    piBEdgeU = (jac(i, j - 1, k + 1) * var.pip(i, j - 1, k) +
                                jac(i, j - 1, k) * var.pip(i, j - 1, k + 1)) /
                               (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    piFEdgeD = (jac(i, j + 1, k - 1) * var.pip(i, j + 1, k) +
                                jac(i, j + 1, k) * var.pip(i, j + 1, k - 1)) /
                               (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
                    piBEdgeD = (jac(i, j - 1, k - 1) * var.pip(i, j - 1, k) +
                                jac(i, j - 1, k) * var.pip(i, j - 1, k - 1)) /
                               (jac(i, j - 1, k) + jac(i, j - 1, k - 1))

                    # Compute pressure gradients.
                    piGradZEdgeU = kappainv * mainv2 * pEdgeU / rhow *
                                   (0.5 * met13EdgeU * (piREdgeU - piLEdgeU) / dx +
                                    0.5 * met23EdgeU * (piFEdgeU - piBEdgeU) / dy +
                                    met33EdgeU *
                                    (var.pip(i, j, k + 1) - var.pip(i, j, k)) / dz)
                    piGradZEdgeD = kappainv * mainv2 * pEdgeD / rhowm *
                                   (0.5 * met13EdgeD * (piREdgeD - piLEdgeD) / dx +
                                    0.5 * met23EdgeD * (piFEdgeD - piBEdgeD) / dy +
                                    met33EdgeD *
                                    (var.pip(i, j, k) - var.pip(i, j, k - 1)) / dz)

                    if k == 1
                        piGradZEdgeD = 0.0
                    elseif k == nz
                        piGradZEdgeU = 0.0
                    end

                    piGrad = 0.5 * (piGradZEdgeU + piGradZEdgeD)
                    facw = 1.0
                    #  if spongeLayer
                    #      facw += dt * kr_sp_w_tfc[i, j, k]
                    #  end

                    buoy = -g_ndim * var.rhop[i, j, k] / rho
                    buoy = 1.0 /
                           (facw + rhostrattfc[i, j, k] / rho * bvsstrattfc[i, j, k] * dt^2) *
                           (-rhostrattfc[i, j, k] / rho *
                            bvsstrattfc[i, j, k] *
                            dt *
                            jac[i, j, k] *
                            (wvrt - dt * piGrad) +
                            facw * buoy +
                            rhostrattfc[i, j, k] / rho *
                            bvsstrattfc[i, j, k] *
                            dt *
                            jac[i, j, k] *
                            facw *
                            0.5 *
                            (met[i, j, k, 1, 3] *
                             (var.u[i, j, k] + var.u[i-1, j, k]) +
                             met[i, j, k, 2, 3] * (var.v[i, j, k] + var.v[i, j-1, k])))

                    var.rhop[i, j, k] = -buoy * rho / g_ndim
                end
                kr_sp_tfc ./= facray
                kr_sp_w_tfc ./= facray
            elseif int_mod == "expl"
                for k in 1:nz, j in 1:ny, i in 1:nx
                    rhop = var.rhop[i, j, k]
                    rho = var.rho[i, j, k] + rhostrattfc[i, j, k]
                    wvrt = 0.5 * (vertWind(i, j, k, model) + vertWind(i, j, k - 1, model))
                    buoy = -g_ndim * rhop / rho
                    buoy -= dt * rhostrattfc[i, j, k] / rho * bvsstrattfc[i, j, k] * wvrt
                    var.rhop[i, j, k] = -buoy * rho / g_ndim
                end
            else
                error("int_mod unknown")
            end
        else
            error("upd_mod unknown")
        end
    else
        error("upd_var unknown")
    end
    end # Timer
end

function momentumPredictor!(model, ode, dt, mmp_mod, int_mod, m, facray)
    @trixi_timeit timer() "Momentum predictor" begin
    #! format: noindent
    # HELP
    (; nx, ny, nz) = model.domain
    (; dx, dy, dz, met, jac) = model.grid
    (; alphark, betark) = ode
    (; kappainv, mainv2, g_ndim) = model.constants
    var = model.variables.prognostic_fields
    flux = model.fluxes
    (; kr_sp_tfc, kr_sp_w_tfc) = model.atmosphere.sponge
    (; rhostrattfc, bvsstrattfc, pstrattfc) = model.atmosphere
    (fluxdiffu, fluxdiffv) = (OffsetArray(zeros((2, 2)), 0:1, 0:1) for _ in 1:2)

    dMom = model.variables.tendencies.dmom
    rhoOld = model.variables.history.rho
    vOld = model.variables.history.v
    uOld = model.variables.history.u
    rhopOld = model.variables.history.rhop

    # TODO
    usave = copy(var.u)

    if m == 1
        dMom .= 0
    end

    if mmp_mod == "rhs"
        if int_mod == "expl"
            #   spongeLayer_s = spongeLayer
            #   spongeLayer = false
        elseif int_mod == "impl"
            kr_sp_tfc .= kr_sp_tfc .* facray
            kr_sp_w_tfc .= kr_sp_w_tfc .* facray
        end
    end

    ################ u*

    i0 = 0
    i1 = nx

    zBoundary = "solid_wall"

    if mmp_mod == "lhs"
        for k in 1:nz
            for j in 1:ny
                for i in i0:i1
                    fL = flux.u[i-1, j, k, 1]
                    fR = flux.u[i, j, k, 1]
                    gB = flux.u[i, j-1, k, 2]
                    gF = flux.u[i, j, k, 2]
                    hD = flux.u[i, j, k-1, 3]
                    hU = flux.u[i, j, k, 3]

                    fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

                    jacEdgeR = 0.5f0 * (jac(i, j, k) + jac(i + 1, j, k))
                    fluxDiff = fluxDiff / jacEdgeR

                    #TODO: Add Coriolis
                    volForce = 0.0
                    if mmp_mod == "lhs"
                        uOld[i, j, k] = var.u[i, j, k]
                        vC = 0.5 * (var.v[i, j, k] + var.v[i, j-1, k])
                        vR = 0.5 * (var.v[i+1, j, k] + var.v[i+1, j-1, k])
                        volForce = volForce +
                                   0.5 *
                                   0.0 * # 0.0 = f_cor_nd[j] <-TODO:
                                   ((rhoOld[i, j, k] + rhostrattfc[i, j, k]) * vC +
                                    (rhoOld[i+1, j, k] + rhostrattfc[i+1, j, k]) * vR)
                        F = -fluxDiff + volForce
                    end

                    rhoM_1 = 0.5 * (rhoOld[i, j, k] + rhoOld[i+1, j, k])
                    rhoM = 0.5 * (var.rho[i, j, k] + var.rho[i+1, j, k])

                    rhoStratEdgeR = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i+1, j, k])
                    rhoM_1 = rhoM_1 + rhoStratEdgeR
                    rhoM = rhoM + rhoStratEdgeR

                    uM_1 = var.u[i, j, k]
                    momM_1 = rhoM_1 * uM_1

                    dMom[i, j, k, 1] = dt * F + alphark[m] * dMom[i, j, k, 1]

                    momM = momM_1 + betark[m] * dMom[i, j, k, 1]

                    uAst = momM / rhoM
                    var.u[i, j, k] = uAst
                end
            end
        end

    elseif mmp_mod == "rhs"
        if int_mod == "expl"
            for k in 1:nz
                for j in 1:ny
                    for i in i0:i1

                        # Compute rhou
                        rhou = 0.5 * (var.rho[i, j, k] + var.rho[i+1, j, k])
                        rhoStratEdgeR = 0.5 *
                                        (rhostrattfc[i, j, k] + rhostrattfc[i+1, j, k])
                        rhou += rhoStratEdgeR

                        # Compute pi values
                        piR = var.pip[i+1, j, k]
                        piL = var.pip[i, j, k]

                        # Compute values at cell edges
                        pEdgeR = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i+1, j, k])
                        met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))

                        # Compute pressure gradient component
                        if k == 1
                            piUUEdgeR = 0.5 * (var.pip[i, j, k+2] +
                                               var.pip[i+1, j, k+2])
                            piUEdgeR = 0.5 *
                                       (var.pip[i, j, k+1] +
                                        var.pip[i+1, j, k+1])
                            piEdgeR = 0.5 *
                                      (var.pip[i, j, k] + var.pip[i+1, j, k])
                            piGrad = kappainv * mainv2 * pEdgeR / rhou *
                                     ((piR - piL) / dx +
                                      met13EdgeR *
                                      (-piUUEdgeR + 4.0 * piUEdgeR - 3.0 * piEdgeR) *
                                      0.5 / dz)

                        elseif k == nz
                            piDDEdgeR = 0.5 * (var.pip[i, j, k-2] +
                                               var.pip[i+1, j, k-2])
                            piDEdgeR = 0.5 *
                                       (var.pip[i, j, k-1] +
                                        var.pip[i+1, j, k-1])
                            piEdgeR = 0.5 *
                                      (var.pip[i, j, k] + var.pip[i+1, j, k])
                            piGrad = kappainv * mainv2 * pEdgeR / rhou *
                                     ((piR - piL) / dx +
                                      met13EdgeR *
                                      (piDDEdgeR - 4.0 * piDEdgeR + 3.0 * piEdgeR) *
                                      0.5 / dz)

                        else
                            piUEdgeR = 0.5 *
                                       (var.pip[i, j, k+1] +
                                        var.pip[i+1, j, k+1])
                            piDEdgeR = 0.5 *
                                       (var.pip[i, j, k-1] +
                                        var.pip[i+1, j, k-1])
                            piGrad = kappainv * mainv2 * pEdgeR / rhou *
                                     ((piR - piL) / dx +
                                      met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                        end

                        volfcx = 0.0

                        # Compute ustar
                        uhorx = var.u[i, j, k]

                        # Coriolis force is integrated on LHS
                        uAst = uhorx + dt * (-piGrad + volfcx / rhou)

                        # if spongeLayer && sponge_uv
                        #     uAst -=
                        #         dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i+1, j, k]) * uhorx
                        # end

                        usave[i, j, k] = uAst
                    end
                end
            end

        elseif int_mod == "impl"
            for k in 1:nz
                for j in 1:ny
                    for i in i0:i1
                        rhou = 0.5 * (var.rho[i, j, k] + var.rho[i+1, j, k])

                        rhov0m = 0.5 * (var.rho[i, j, k] + var.rho[i, j-1, k])
                        rhov00 = 0.5 * (var.rho[i, j+1, k] + var.rho[i, j, k])
                        rhov1m = 0.5 * (var.rho[i+1, j, k] + var.rho[i+1, j-1, k])
                        rhov10 = 0.5 * (var.rho[i+1, j+1, k] + var.rho[i+1, j, k])

                        rhoStratEdgeR = 0.5 *
                                        (rhostrattfc[i, j, k] + rhostrattfc[i+1, j, k])
                        rhou += rhoStratEdgeR

                        piR = var.pip[i+1, j, k]
                        piL = var.pip[i, j, k]

                        pEdgeR = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i+1, j, k])
                        met13EdgeR = 0.5 * (met[i, j, k, 1, 3] + met[i+1, j, k, 1, 3])

                        if k == 1 && zBoundary == "solid_wall"
                            piUUEdgeR = 0.5 * (var.pip[i, j, k+2] +
                                               var.pip[i+1, j, k+2])
                            piUEdgeR = 0.5 *
                                       (var.pip[i, j, k+1] +
                                        var.pip[i+1, j, k+1])
                            piEdgeR = 0.5 *
                                      (var.pip[i, j, k] + var.pip[i+1, j, k])
                            piGradX = kappainv * mainv2 * pEdgeR / rhou *
                                      ((piR - piL) / dx +
                                       met13EdgeR *
                                       (-piUUEdgeR + 4.0 * piUEdgeR - 3.0 * piEdgeR) *
                                       0.5 / dz)
                        elseif k == nz && zBoundary == "solid_wall"
                            piDDEdgeR = 0.5 * (var.pip[i, j, k-2] +
                                               var.pip[i+1, j, k-2])
                            piDEdgeR = 0.5 *
                                       (var.pip[i, j, k-1] +
                                        var.pip[i+1, j, k-1])
                            piEdgeR = 0.5 *
                                      (var.pip[i, j, k] + var.pip[i+1, j, k])
                            piGradX = kappainv * mainv2 * pEdgeR / rhou *
                                      ((piR - piL) / dx +
                                       met13EdgeR *
                                       (piDDEdgeR - 4.0 * piDEdgeR + 3.0 * piEdgeR) *
                                       0.5 / dz)
                        else
                            piUEdgeR = 0.5 *
                                       (var.pip[i, j, k+1] +
                                        var.pip[i+1, j, k+1])
                            piDEdgeR = 0.5 *
                                       (var.pip[i, j, k-1] +
                                        var.pip[i+1, j, k-1])
                            piGradX = kappainv * mainv2 * pEdgeR / rhou *
                                      ((piR - piL) / dx +
                                       met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                        end

                        volfcx = 0.0
                        volfcy = 0.0

                        uhorx = var.u[i, j, k]

                        vhory = 0.25 * (var.v[i, j-1, k] +
                                        var.v[i, j, k] +
                                        var.v[i+1, j-1, k] +
                                        var.v[i+1, j, k])

                        facu = 1.0

                        #  if spongeLayer && sponge_uv
                        #      facu += dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k])
                        #  end

                        facv = facu

                        uAst = 1.0 / facu * (uhorx + dt * (-piGradX + volfcx / rhou))
                        usave[i, j, k] = uAst
                    end
                end
            end
        end
    end

    #####
    # V correction
    j0 = 0
    j1 = ny
    if mmp_mod == "lhs"
        if int_mod != "expl"
            error("ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs")
        end

        for k in 1:nz, j in j0:j1, i in 1:nx
            fR = flux.v[i, j, k, 1]
            fL = flux.v[i-1, j, k, 1]
            gF = flux.v[i, j, k, 2]
            gB = flux.v[i, j-1, k, 2]
            hU = flux.v[i, j, k, 3]
            hD = flux.v[i, j, k-1, 3]
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            jacEdgeF = 0.5 * (jac[i, j, k] + jac[i, j+1, k])
            fluxDiff /= jacEdgeF
            volForce = 0.0

            if mmp_mod == "lhs"
                vOld[i, j, k] = var.v[i, j, k]
                uC = 0.5 * (uOld[i, j, k] + uOld[i-1, j, k])
                uF = 0.5 * (uOld[i, j+1, k] + uOld[i-1, j+1, k])
                #   volForce -= 0.5 * (f_cor_nd[j] * (rhoOld[i, j, k] + rhoStrat[i, j, k]) * uC +
                #                     f_cor_nd[j + 1] * (rhoOld[i, j + 1, k] + rhoStrat[i, j + 1, k]) * uF)
            end

            if mmp_mod == "lhs"
                F = -fluxDiff + volForce
            else
                error("ERROR: wrong mmp_mod")
            end

            # if model == "pseudo_incompressible"
            rhoM_1 = 0.5 * (rhoOld[i, j, k] + rhoOld[i, j+1, k])
            rhoM = 0.5 * (var.rho[i, j, k] + var.rho[i, j+1, k])
            rhoStratEdgeF = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j+1, k])
            rhoM_1 += rhoStratEdgeF
            rhoM += rhoStratEdgeF
            #  else
            #  error("momentumPredictor: unknown case model.")
            #   end

            vM_1 = var.v[i, j, k]
            momM_1 = rhoM_1 * vM_1
            dMom[i, j, k, 2] = dt * F + alphark[m] * dMom[i, j, k, 2]
            momM = momM_1 + betark[m] * dMom[i, j, k, 2]
            vAst = momM / rhoM
            var.v[i, j, k] = vAst
        end
    elseif mmp_mod == "rhs"
        if int_mod == "expl"
            for k in 1:nz, j in j0:j1, i in 1:nx
                rhov = 0.5 * (var.rho[i, j, k] + var.rho[i, j+1, k])
                rhoStratEdgeF = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j+1, k])
                rhov += rhoStratEdgeF

                piF = var.pip[i, j+1, k]
                piB = var.pip[i, j, k]
                pEdgeF = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j+1, k])
                met23EdgeF = 0.5 * (met[i, j, k, 2, 3] + met[i, j+1, k, 2, 3])

                if k == 1 && zBoundary == "solid_wall"
                    piUUEdgeF = 0.5 *
                                (var.pip[i, j, k+2] + var.pip[i, j+1, k+2])
                    piUEdgeF = 0.5 *
                               (var.pip[i, j, k+1] + var.pip[i, j+1, k+1])
                    piEdgeF = 0.5 * (var.pip[i, j, k] + var.pip[i, j+1, k])
                    piGrad = kappainv * mainv2 * pEdgeF / rhov * ((piF - piB) / dy +
                                                                  met23EdgeF *
                                                                  (-piUUEdgeF + 4.0 * piUEdgeF - 3.0 * piEdgeF) *
                                                                  0.5 / dz)
                elseif k == nz && zBoundary == "solid_wall"
                    piDDEdgeF = 0.5 *
                                (var.pip[i, j, k-2] + var.pip[i, j+1, k-2])
                    piDEdgeF = 0.5 *
                               (var.pip[i, j, k-1] + var.pip[i, j+1, k-1])
                    piEdgeF = 0.5 * (var.pip[i, j, k] + var.pip[i, j+1, k])
                    piGrad = kappainv * mainv2 * pEdgeF / rhov * ((piF - piB) / dy +
                                                                  met23EdgeF *
                                                                  (piDDEdgeF - 4.0 * piDEdgeF + 3.0 * piEdgeF) *
                                                                  0.5 / dz)
                else
                    piUEdgeF = 0.5 *
                               (var.pip[i, j, k+1] + var.pip[i, j+1, k+1])
                    piDEdgeF = 0.5 *
                               (var.pip[i, j, k-1] + var.pip[i, j+1, k-1])
                    piGrad = kappainv * mainv2 * pEdgeF / rhov *
                             ((piF - piB) / dy +
                              met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                end

                vAst = var.v[i, j, k] + dt * (-piGrad)
                ## TODO:
                #  if spongeLayer && sponge_uv
                #     vAst -= dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) * var.v[i, j, k]
                # end
                var.v[i, j, k] = vAst
            end
        elseif int_mod == "impl"
            for k in 1:nz, j in j0:j1, i in 1:nx
                rhoum0 = 0.5 * (var.rho[i, j, k] + var.rho[i-1, j, k])
                rhou00 = 0.5 * (var.rho[i+1, j, k] + var.rho[i, j, k])
                rhoum1 = 0.5 * (var.rho[i, j+1, k] + var.rho[i-1, j+1, k])
                rhou01 = 0.5 * (var.rho[i+1, j+1, k] + var.rho[i, j+1, k])

                rhov = 0.5 * (var.rho[i, j, k] + var.rho[i, j+1, k])

                rhoStratEdgeF = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i, j+1, k])
                rhov += rhoStratEdgeF

                piF = var.pip[i, j+1, k]
                piB = var.pip[i, j, k]

                pEdgeF = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j+1, k])
                met23EdgeF = 0.5 * (met[i, j, k, 2, 3] + met[i, j+1, k, 2, 3])

                if k == 1 && zBoundary == "solid_wall"
                    piUUEdgeF = 0.5 *
                                (var.pip[i, j, k+2] + var.pip[i, j+1, k+2])
                    piUEdgeF = 0.5 *
                               (var.pip[i, j, k+1] + var.pip[i, j+1, k+1])
                    piEdgeF = 0.5 * (var.pip[i, j, k] + var.pip[i, j+1, k])
                    piGradY = kappainv * mainv2 * pEdgeF / rhov * ((piF - piB) / dy +
                                                                   met23EdgeF *
                                                                   (-piUUEdgeF + 4.0 * piUEdgeF - 3.0 * piEdgeF) *
                                                                   0.5 / dz)
                elseif k == nz && zBoundary == "solid_wall"
                    piDDEdgeF = 0.5 *
                                (var.pip[i, j, k-2] + var.pip[i, j+1, k-2])
                    piDEdgeF = 0.5 *
                               (var.pip[i, j, k-1] + var.pip[i, j+1, k-1])
                    piEdgeF = 0.5 * (var.pip[i, j, k] + var.pip[i, j+1, k])
                    piGradY = kappainv * mainv2 * pEdgeF / rhov * ((piF - piB) / dy +
                                                                   met23EdgeF *
                                                                   (piDDEdgeF - 4.0 * piDEdgeF + 3.0 * piEdgeF) *
                                                                   0.5 / dz)
                else
                    piUEdgeF = 0.5 *
                               (var.pip[i, j, k+1] + var.pip[i, j+1, k+1])
                    piDEdgeF = 0.5 *
                               (var.pip[i, j, k-1] + var.pip[i, j+1, k-1])
                    piGradY = kappainv * mainv2 * pEdgeF / rhov *
                              ((piF - piB) / dy +
                               met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                end

                volfcx = 0.0
                volfcy = 0.0

                uhorx = 0.25 * (var.u[i-1, j, k] +
                                var.u[i-1, j+1, k] +
                                var.u[i, j, k] +
                                var.u[i, j+1, k])
                vhory = var.v[i, j, k]

                facv = 1.0
                #if spongeLayer && sponge_uv
                #   facv += dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k])
                #end
                #@assert false

                facu = facv

                vAst = 1.0 / facv * (vhory + dt * (-piGradY + volfcy / rhov))

                var.v[i, j, k] = vAst
            end
        else
            error("ERROR: unknown int_mod")
        end
        var.u .= usave
    else
        error("ERROR: unknown mmp_mod")
    end
    k0 = 1
    k1 = nz - 1
    ################## W Correction
    if mmp_mod == "lhs"
        for k in k0:k1, j in 1:ny, i in 1:nx

            # Compute vertical momentum flux divergence.
            fR = flux.w[i, j, k, 1]
            fL = flux.w[i-1, j, k, 1]
            gF = flux.w[i, j, k, 2]
            gB = flux.w[i, j-1, k, 2]
            hU = flux.w[i, j, k, 3]
            hD = flux.w[i, j, k-1, 3]

            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
            # Adjust Cartesian vertical momentum flux divergence.
            jacEdgeU = 2.0 * jac[i, j, k] * jac[i, j, k+1] /
                       (jac[i, j, k] + jac[i, j, k+1])
            fluxDiff = fluxDiff / jacEdgeU

            # Compute zonal momentum flux divergences.
            for ll in 0:1
                for mm in 0:1
                    fR = flux.u[i-ll, j, k+mm, 1]
                    fL = flux.u[i-1-ll, j, k+mm, 1]
                    gF = flux.u[i-ll, j, k+mm, 2]
                    gB = flux.u[i-ll, j-1, k+mm, 2]
                    hU = flux.u[i-ll, j, k+mm, 3]
                    hD = flux.u[i-ll, j, k-1+mm, 3]
                    fluxdiffu[ll, mm] = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
                    jacEdgeR = 0.5 *
                               (jac[i-ll, j, k+mm] + jac[i+1-ll, j, k+mm])
                    fluxdiffu[ll, mm] = fluxdiffu[ll, mm] / jacEdgeR
                end
            end

            # Compute meridional momentum flux divergences.
            for ll in 0:1
                for mm in 0:1
                    fR = flux.v[i, j-ll, k+mm, 1]
                    fL = flux.v[i-1, j-ll, k+mm, 1]
                    gF = flux.v[i, j-ll, k+mm, 2]
                    gB = flux.v[i, j-1-ll, k+mm, 2]
                    hU = flux.v[i, j-ll, k+mm, 3]
                    hD = flux.v[i, j-ll, k-1+mm, 3]
                    fluxdiffv[ll, mm] = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
                    jacEdgeF = 0.5 *
                               (jac[i, j-ll, k+mm] + jac[i, j+1-ll, k+mm])
                    fluxdiffv[ll, mm] = fluxdiffv[ll, mm] / jacEdgeF
                end
            end

            # Compute transformed vertical momentum flux divergence.
            fluxDiff = trafo(i,
                j,
                k,
                fluxdiffu[0, 0],
                fluxdiffu[0, 1],
                fluxdiffu[1, 0],
                fluxdiffu[1, 1],
                fluxdiffv[0, 0],
                fluxdiffv[0, 1],
                fluxdiffv[1, 0],
                fluxdiffv[1, 1],
                fluxDiff,
                "tfc",
                model)

            # Explicit integration of Coriolis force in TFC.
            if mmp_mod == "lhs"
                #    vC = 0.5 * (vOldTFC[i, j, k] + vOldTFC[i, j - 1, k])
                #    vU = 0.5 * (vOldTFC[i, j, k + 1] + vOldTFC[i, j - 1, k + 1])
                #    uC = 0.5 * (uOldTFC[i, j, k] + uOldTFC[i - 1, j, k])
                #    uU = 0.5 * (uOldTFC[i, j, k + 1] + uOldTFC[i - 1, j, k + 1])

                #    volForce = volForce + f_cor_nd[j] * (jac[i, j, k + 1] * met[i, j, k, 1, 3] * (rhoOld[i, j, k] + rhoStratTFC[i, j, k]) * vC + jac[i, j, k] * met[i, j, k + 1, 1, 3] * (rhoOld[i, j, k + 1] + rhoStratTFC[i, j, k + 1]) * vU) / (jac[i, j, k] + jac[i, j, k + 1]) - f_cor_nd[j] * (jac[i, j, k + 1] * met[i, j, k, 2, 3] * (rhoOld[i, j, k] + rhoStratTFC[i, j, k]) * uC + jac[i, j, k] * met[i, j, k + 1, 2, 3] * (rhoOld[i, j, k + 1] + rhoStratTFC[i, j, k + 1]) * uU) / (jac[i, j, k] + jac[i, j, k + 1])
            end

            #--------------------
            #   F(phi) = RHS
            #--------------------
            volForce = 0.0
            if mmp_mod == "lhs"
                F = -fluxDiff + volForce
            else
                # error("ERROR: wrong mmp_mod")
            end

            # interpolated densities
            #if model == "pseudo_incompressible"
            rhoM_1 = (jac[i, j, k+1] * rhoOld[i, j, k] +
                      jac[i, j, k] * rhoOld[i, j, k+1]) /
                     (jac[i, j, k] + jac[i, j, k+1])
            rhoM = (jac[i, j, k+1] * var.rho[i, j, k] +
                    jac[i, j, k] * var.rho[i, j, k+1]) /
                   (jac[i, j, k] + jac[i, j, k+1])

            rhoStratEdgeU = (jac[i, j, k+1] * rhostrattfc[i, j, k] +
                             jac[i, j, k] * rhostrattfc[i, j, k+1]) /
                            (jac[i, j, k] + jac[i, j, k+1])
            rhoM_1 = rhoM_1 + rhoStratEdgeU
            rhoM = rhoM + rhoStratEdgeU
            #else
            # error("momentumPredictor: unknown case model.")
            #end

            # velocity and momentum at t(m-1)
            wM_1 = var.w[i, j, k]
            momM_1 = rhoM_1 * wM_1

            # q(m-1) -> q(m)
            dMom[i, j, k, 3] = dt * F + alphark[m] * dMom[i, j, k, 3]

            # rhoW(m-1) -> rhoW(m)
            momM = momM_1 + betark[m] * dMom[i, j, k, 3]

            # calc w(m,*)
            wAst = momM / rhoM

            # wAst -> var
            var.w[i, j, k] = wAst
        end
    elseif mmp_mod == "rhs"
        if int_mod == "expl"
            for k in k0:k1, j in 1:ny, i in 1:nx
                rho000 = var.rho[i, j, k]
                rho001 = var.rho[i, j, k+1]

                rhow = (jac[i, j, k+1] * var.rho[i, j, k] +
                        jac[i, j, k] * var.rho[i, j, k+1]) /
                       (jac[i, j, k] + jac[i, j, k+1])

                rhoStratEdgeU = (jac[i, j, k+1] * rhostrattfc[i, j, k] +
                                 jac[i, j, k] * rhostrattfc[i, j, k+1]) /
                                (jac[i, j, k] + jac[i, j, k+1])
                rho000 = rho000 + rhostrattfc[i, j, k]
                rho001 = rho001 + rhostrattfc[i, j, k+1]
                rhow = rhow + rhoStratEdgeU

                piU = var.pip[i, j, k+1]
                piD = var.pip[i, j, k]

                # Compute values at cell edges.
                pEdgeU = (jac[i, j, k+1] * pstrattfc[i, j, k] +
                          jac[i, j, k] * pstrattfc[i, j, k+1]) /
                         (jac[i, j, k] + jac[i, j, k+1])
                met13EdgeU = (jac[i, j, k+1] * met[i, j, k, 1, 3] +
                              jac[i, j, k] * met[i, j, k+1, 1, 3]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                met23EdgeU = (jac[i, j, k+1] * met[i, j, k, 2, 3] +
                              jac[i, j, k] * met[i, j, k+1, 2, 3]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                met33EdgeU = (jac[i, j, k+1] * met[i, j, k, 3, 3] +
                              jac[i, j, k] * met[i, j, k+1, 3, 3]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                piREdgeU = (jac[i+1, j, k+1] * var.pip[i+1, j, k] +
                            jac[i+1, j, k] * var.pip[i+1, j, k+1]) /
                           (jac[i+1, j, k] + jac[i+1, j, k+1])
                piLEdgeU = (jac[i-1, j, k+1] * var.pip[i-1, j, k] +
                            jac[i-1, j, k] * var.pip[i-1, j, k+1]) /
                           (jac[i-1, j, k] + jac[i-1, j, k+1])
                piFEdgeU = (jac[i, j+1, k+1] * var.pip[i, j+1, k] +
                            jac[i, j+1, k] * var.pip[i, j+1, k+1]) /
                           (jac[i, j+1, k] + jac[i, j+1, k+1])
                piBEdgeU = (jac[i, j-1, k+1] * var.pip[i, j-1, k] +
                            jac[i, j-1, k] * var.pip[i, j-1, k+1]) /
                           (jac[i, j-1, k] + jac[i, j-1, k+1])

                # Compute pressure gradient component.
                piGrad = kappainv * mainv2 * pEdgeU / rhow *
                         (met13EdgeU * (piREdgeU - piLEdgeU) * 0.5 / dx +
                          met23EdgeU * (piFEdgeU - piBEdgeU) * 0.5 / dy +
                          met33EdgeU * (var.pip[i, j, k+1] - var.pip[i, j, k]) /
                          dz)

                volfcz = 0.0

                # wstar
                wvert = var.w[i, j, k]

                buoy = -g_ndim *
                       (jac[i, j, k+1] * rhopOld[i, j, k] / rho000 / jac[i, j, k] +
                        jac[i, j, k] * rhopOld[i, j, k+1] / rho001 / jac[i, j, k+1]) /
                       (jac[i, j, k] + jac[i, j, k+1])

                wAst = wvert + dt * (buoy - piGrad + volfcz / rhow)

                #if spongeLayer
                #    wAst = wAst - dt * (jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] + jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]) / (jac[i, j, k] + jac[i, j, k + 1]) * wvert
                #end

                var.w[i, j, k] = wAst
            end
        elseif int_mod == "impl"
            for k in k0:k1, j in 1:ny, i in 1:nx
                rho000 = var.rho[i, j, k]
                rho001 = var.rho[i, j, k+1]

                rhow = (jac[i, j, k+1] * var.rho[i, j, k] +
                        jac[i, j, k] * var.rho[i, j, k+1]) /
                       (jac[i, j, k] + jac[i, j, k+1])

                rhoStratEdgeU = (jac[i, j, k+1] * rhostrattfc[i, j, k] +
                                 jac[i, j, k] * rhostrattfc[i, j, k+1]) /
                                (jac[i, j, k] + jac[i, j, k+1])
                rho000 = rho000 + rhostrattfc[i, j, k]
                rho001 = rho001 + rhostrattfc[i, j, k+1]
                rhow = rhow + rhoStratEdgeU

                piU = var.pip[i, j, k+1]
                piD = var.pip[i, j, k]

                # Compute values at cell edges.
                pEdgeU = (jac[i, j, k+1] * pstrattfc[i, j, k] +
                          jac[i, j, k] * pstrattfc[i, j, k+1]) /
                         (jac[i, j, k] + jac[i, j, k+1])
                met13EdgeU = (jac[i, j, k+1] * met[i, j, k, 1, 3] +
                              jac[i, j, k] * met[i, j, k+1, 1, 3]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                met23EdgeU = (jac[i, j, k+1] * met[i, j, k, 2, 3] +
                              jac[i, j, k] * met[i, j, k+1, 2, 3]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                met33EdgeU = (jac[i, j, k+1] * met[i, j, k, 3, 3] +
                              jac[i, j, k] * met[i, j, k+1, 3, 3]) /
                             (jac[i, j, k] + jac[i, j, k+1])
                piREdgeU = (jac[i+1, j, k+1] * var.pip[i+1, j, k] +
                            jac[i+1, j, k] * var.pip[i+1, j, k+1]) /
                           (jac[i+1, j, k] + jac[i+1, j, k+1])
                piLEdgeU = (jac[i-1, j, k+1] * var.pip[i-1, j, k] +
                            jac[i-1, j, k] * var.pip[i-1, j, k+1]) /
                           (jac[i-1, j, k] + jac[i-1, j, k+1])
                piFEdgeU = (jac[i, j+1, k+1] * var.pip[i, j+1, k] +
                            jac[i, j+1, k] * var.pip[i, j+1, k+1]) /
                           (jac[i, j+1, k] + jac[i, j+1, k+1])
                piBEdgeU = (jac[i, j-1, k+1] * var.pip[i, j-1, k] +
                            jac[i, j-1, k] * var.pip[i, j-1, k+1]) /
                           (jac[i, j-1, k] + jac[i, j-1, k+1])

                # Compute pressure gradient component.
                piGrad = kappainv * mainv2 * pEdgeU / rhow *
                         (met13EdgeU * (piREdgeU - piLEdgeU) * 0.5 / dx +
                          met23EdgeU * (piFEdgeU - piBEdgeU) * 0.5 / dy +
                          met33EdgeU * (var.pip[i, j, k+1] - var.pip[i, j, k]) /
                          dz)

                volfcz = 0.0

                # wstar
                wvert = var.w[i, j, k]

                # squared Brunt-Vaisala frequency averaged to half levels
                bvsstw = (jac[i, j, k+1] * bvsstrattfc[i, j, k] +
                          jac[i, j, k] * bvsstrattfc[i, j, k+1]) /
                         (jac[i, j, k] + jac[i, j, k+1])

                facw = 1.0

                if model.parameters.boundaries.spongelayer == true
                    facw = facw +
                           dt * (jac[i, j, k+1] * kr_sp_w_tfc[i, j, k] +
                                 jac[i, j, k] * kr_sp_w_tfc[i, j, k+1]) /
                           (jac[i, j, k] + jac[i, j, k+1])
                end

                # Buoyancy is predicted after momentum in implicit steps.
                buoy = -g_ndim *
                       (jac[i, j, k+1] * var.rhop[i, j, k] / rho000 / jac[i, j, k] +
                        jac[i, j, k] * var.rhop[i, j, k+1] / rho001 /
                        jac[i, j, k+1]) /
                       (jac[i, j, k] + jac[i, j, k+1])

                uC = 0.5 * (var.u[i, j, k] + var.u[i-1, j, k])
                uU = 0.5 * (var.u[i, j, k+1] + var.u[i-1, j, k+1])
                vC = 0.5 * (var.v[i, j, k] + var.v[i, j-1, k])
                vU = 0.5 * (var.v[i, j, k+1] + var.v[i, j-1, k+1])

                wAst = 1.0 / (facw + rhoStratEdgeU / rhow * bvsstw * dt^2.0) *
                       (wvert - dt * piGrad +
                        dt * buoy +
                        dt * volfcz / rhow +
                        rhoStratEdgeU / rhow *
                        bvsstw *
                        dt^2.0 *
                        (jac[i, j, k+1] *
                         (met[i, j, k, 1, 3] * uC + met[i, j, k, 2, 3] * vC) +
                         jac[i, j, k] *
                         (met[i, j, k+1, 1, 3] * uU + met[i, j, k+1, 2, 3] * vU)) /
                        (jac[i, j, k] + jac[i, j, k+1]))

                var.w[i, j, k] = wAst
            end
        end
    end

    if mmp_mod == "rhs"
        if int_mod == "expl"
            #  spongeLayer = spongeLayer_s
        elseif int_mod == "impl"
            kr_sp_tfc .= kr_sp_tfc ./ facray
            kr_sp_w_tfc .= kr_sp_w_tfc ./ facray
        end
    end
    end # timer
end #END
