function compute_fluxDiff(flux, i, j, k, grid)
    (; dx, dy, dz) = grid
    fL = flux[i-1, j, k, 1]
    fR = flux[i, j, k, 1]
    gB = flux[i, j-1, k, 2]
    gF = flux[i, j, k, 2]
    hD = flux[i, j, k-1, 3]
    hU = flux[i, j, k, 3]

    fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

    return fluxDiff
end

function massUpdate_rho!(semi, ode, upd_mod, m)

    (; grid, cache, equations) = semi
    (; var, flux, rhoStrat, bvsStrat, jac) = cache
    (; nx, ny, nz, lz) = grid
    (; dt, alphaRK, betaRK, dRho, dRhop) = ode
    (; g_ndim) = equations

    # I think m is the temporal level or stage
    if m == 1
        dRho .= 0.0
    end

    for k = 1:nz
        for j = 1:ny
            for i = 1:nx

                fluxDiff = compute_fluxDiff(flux.rho, i, j, k, grid)

                fluxDiff = fluxDiff / jac(i, j, k)

                F = -fluxDiff

                dRho[i, j, k] = dt * F + alphaRK[m] * dRho[i, j, k]

                var.rho[i, j, k] = var.rho[i, j, k] + betaRK[m] * dRho[i, j, k]
            end
        end
    end

end


function massUpdate_rhop_impl!(semi, ode, upd_mod, m)

    (; grid, cache, equations) = semi
    (; var, flux, rhoStrat, bvsStrat, jac) = cache
    (; nx, ny, nz, lz) = grid
    (; dt, alphaRK, betaRK, dRho, dRhop) = ode
    (; g_ndim) = equations

    # I think m is the temporal level or stage
    if m == 1
        dRhop .= 0.0
    end
    
    zBoundary = "solid_wall"

    kr_sp_tfc *= facray
    kr_sp_w_tfc *= facray
    
    for k in 1:nz, j in 1:ny, i in 1:nx
        rho = var.rho[i, j, k]
        
        rhow = (jac(i, j, k + 1) * var.rho[i, j, k] + jac(i, j, k) * var.rho[i, j, k + 1]) / (jac(i, j, k) + jac(i, j, k + 1))
        rhowm = (jac(i, j, k - 1) * var.rho[i, j, k] + jac(i, j, k) * var.rho[i, j, k - 1]) / (jac(i, j, k) + jac(i, j, k - 1))
        
        rho += rhoStratTFC[i, j, k]
        rhow += (jac(i, j, k + 1) * rhoStratTFC[i, j, k] + jac(i, j, k) * rhoStratTFC[i, j, k + 1]) / (jac(i, j, k) + jac(i, j, k + 1))
        rhowm += (jac(i, j, k - 1) * rhoStratTFC[i, j, k] + jac(i, j, k) * rhoStratTFC[i, j, k - 1]) / (jac(i, j, k) + jac(i, j, k - 1))
        
        wvrt = 0.5 * (wOldTFC[i, j, k] + wOldTFC[i, j, k - 1])
        
        pEdgeU = (jac(i, j, k + 1) * pStratTFC[i, j, k] + jac(i, j, k) * pStratTFC[i, j, k + 1]) / (jac(i, j, k) + jac(i, j, k + 1))
        pEdgeD = (jac(i, j, k - 1) * pStratTFC[i, j, k] + jac(i, j, k) * pStratTFC[i, j, k - 1]) / (jac(i, j, k) + jac(i, j, k - 1))
        
        met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
        met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
        met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
        met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) + jac(i, j, k) * met(i, j, k - 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
        met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) + jac(i, j, k) * met(i, j, k - 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
        met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) + jac(i, j, k) * met(i, j, k - 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
        
        piGradZEdgeU = kappaInv * MaInv2 * pEdgeU / rhow * (0.5 * met13EdgeU * (var.pi[i + 1, j, k] - var.pi[i - 1, j, k]) / dx +
            0.5 * met23EdgeU * (var.pi[i, j + 1, k] - var.pi[i, j - 1, k]) / dy + met33EdgeU * (var.pi[i, j, k + 1] - var.pi[i, j, k]) / dz)
        
        piGradZEdgeD = kappaInv * MaInv2 * pEdgeD / rhowm * (0.5 * met13EdgeD * (var.pi[i + 1, j, k] - var.pi[i - 1, j, k]) / dx +
            0.5 * met23EdgeD * (var.pi[i, j + 1, k] - var.pi[i, j - 1, k]) / dy + met33EdgeD * (var.pi[i, j, k] - var.pi[i, j, k - 1]) / dz)
        
        if k == 1 && zBoundary == "solid_wall"
            piGradZEdgeD = 0.0
        elseif k == nz && zBoundary == "solid_wall"
            piGradZEdgeU = 0.0
        end
        
        piGrad = 0.5 * (piGradZEdgeU + piGradZEdgeD)
        
        facw = 1.0
        if spongeLayer
            facw += dt * kr_sp_w_tfc[i, j, k]
        end
        
        buoy = -g_ndim * var.rhop[i, j, k] / rho
        buoy = 1.0 / (facw + rhoStratTFC[i, j, k] / rho * bvsStratTFC[i, j, k] * dt^2) *
            (-rhoStratTFC[i, j, k] / rho * bvsStratTFC[i, j, k] * dt * jac(i, j, k) * (wvrt - dt * piGrad) +
            facw * buoy + rhoStratTFC[i, j, k] / rho * bvsStratTFC[i, j, k] * dt * jac(i, j, k) * facw *
            0.5 * (met(i, j, k, 1, 3) * (var.u[i, j, k] + var.u[i - 1, j, k]) +
                   met(i, j, k, 2, 3) * (var.v[i, j, k] + var.v[i, j - 1, k])))
        
        var.rhop[i, j, k] = -buoy * rho / g_ndim
    end
    
    kr_sp_tfc /= facray
    kr_sp_w_tfc /= facray

end




function massUpdate_rhop!(semi, ode, upd_mod, m)

    (; grid, cache, equations) = semi
    (; var, flux, rhoStrat, bvsStrat, jac) = cache
    (; nx, ny, nz, lz) = grid
    (; dt, alphaRK, betaRK, dRho, dRhop) = ode
    (; g_ndim) = equations

    # I think m is the temporal level or stage
    if m == 1
        dRhop .= 0.0
    end

    if upd_mod == "lhs"

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    fluxDiff = compute_fluxDiff(flux.rhop, i, j, k, grid)

                    fluxDiff = fluxDiff / jac(i, j, k)

                    F = -fluxDiff

                    dRhop[i, j, k] = dt * F + alphaRK[m] * dRhop[i, j, k]

                    var.rhop[i, j, k] = var.rhop[i, j, k] + betaRK[m] * dRhop[i, j, k]
                end
            end
        end
    elseif upd_mod == "rhs"
        ## LHS
        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    rhop = var.rhop[i, j, k]
                    rho = var.rho[i, j, k]
                    rho = rho + rhoStrat[i, j, k]

                    wvrt =
                        0.5 *
                        (vertWind(i, j, k, semi) + vertWind(i, j, k - 1, semi))

                    buoy = -g_ndim * rhop / rho
                    buoy = buoy - dt * rhoStrat[i, j, k] / rho * bvsStrat[i, j, k] * wvrt
                    var.rhop[i, j, k] = -buoy * rho / g_ndim
                end
            end
        end
    end
end




function momentumPredictor!(semi, ode, mmp_mod, m)

    momentumPredictor_u!(semi, ode, mmp_mod, m)
    momentumPredictor_v!(semi, ode, mmp_mod, m)
    momentumPredictor_w!(semi, ode, mmp_mod, m)

    # if lhs or rhs reuse "saved" variables
end

function momentumPredictor_u!(semi, ode, mmp_mod, m)
    ## unpack boundaries and set the boundaries for x direction
    # periodic case

    (; grid, cache, met, equations) = semi
    (; var, flux, rhoStrat, bvsStrat, pStrat, jac) = cache
    (; nx, ny, nz, lz, dx, dy, dz) = grid
    (; dt, alphaRK, betaRK, dMom, usave, rhoOld) = ode
    (; kappaInv, MaInv2) = equations

    if m == 1
        dMom[:,:,:,1] .= 0
    end

    i0 = 0
    i1 = nx
    if mmp_mod == "lhs"
        for k = 1:nz
            for j = 1:ny
                for i = i0:i1
                    fluxDiff = compute_fluxDiff(flux.u, i, j, k, grid)

                    jacEdgeR = 0.5f0 * (jac(i, j, k) + jac(i + 1, j, k))
                    fluxDiff = fluxDiff / jacEdgeR

                    #TODO: Add Coriolis
                    volForce = 0.0
                    # if mmp_mmod "lhs"
                   # uOld[i, j, k] = var.u[i, j, k]
                   # vC = 0.5 * (var.v[i, j, k] + var.v[i, j-1, k])
                   # vR = 0.5 * (var.v[i+1, j, k] + var.v[i+1, j-1, k])
                    #volForce =
                   #      volForce +
                   #      0.5 *
                   #      f_cor_nd[j] *
                   #      (
                   #          (rhoOld[i, j, k] + rhoStrat[i, j, k]) *
                  #          (rhoOld[i+1, j, k] + rhoStrat[i+1, j, k])
                  #      ) *
                  #      vR
                    F = -fluxDiff + volForce
                    # end


                    rhoM_1 = 0.5 * (rhoOld[i, j, k] + rhoOld[i+1, j, k])
                    rhoM = 0.5 * (var.rho[i, j, k] + var.rho[i+1, j, k])

                    rhoStratEdgeR = 0.5 * (rhoStrat[i, j, k] + rhoStrat[i+1, j, k])
                    rhoM_1 = rhoM_1 + rhoStratEdgeR
                    rhoM = rhoM + rhoStratEdgeR

                    uM_1 = var.u[i, j, k]
                    momM_1 = rhoM_1 * uM_1

                    dMom[i, j, k, 1] = dt * F + alphaRK[m] * dMom[i, j, k, 1]

                    momM = momM_1 + betaRK[m] * dMom[i, j, k, 1]

                    uAst = momM / rhoM
                    var.u[i, j, k] = uAst
                end
            end
        end
    elseif mmp_mod == "rhs"
        # RHS Explicit

        for k = 1:nz
            for j = 1:ny
                for i = i0:i1

                    # Compute rhou
                    rhou = 0.5 * (var.rho[i, j, k] + var.rho[i+1, j, k])
                    rhoStratEdgeR = 0.5 * (rhoStrat[i, j, k] + rhoStrat[i+1, j, k])
                    rhou += rhoStratEdgeR

                    # Compute pi values
                    piR = var.exner[i+1, j, k]
                    piL = var.exner[i, j, k]

                    # Compute values at cell edges
                    pEdgeR = 0.5 * (pStrat[i, j, k] + pStrat[i+1, j, k])
                    met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i+1, j, k, 1, 3))

                    # Compute pressure gradient component
                    if k == 1
                        piUUEdgeR = 0.5 * (var.exner[i, j, k+2] + var.exner[i+1, j, k+2])
                        piUEdgeR = 0.5 * (var.exner[i, j, k+1] + var.exner[i+1, j, k+1])
                        piEdgeR = 0.5 * (var.exner[i, j, k] + var.exner[i+1, j, k])
                        piGrad =
                            kappaInv * MaInv2 * pEdgeR / rhou * (
                                (piR - piL) / dx +
                                met13EdgeR *
                                (-piUUEdgeR + 4.0 * piUEdgeR - 3.0 * piEdgeR) *
                                0.5 / dz
                            )

                    elseif k == nz
                        piDDEdgeR = 0.5 * (var.exner[i, j, k-2] + var.exner[i+1, j, k-2])
                        piDEdgeR = 0.5 * (var.exner[i, j, k-1] + var.exner[i+1, j, k-1])
                        piEdgeR = 0.5 * (var.exner[i, j, k] + var.exner[i+1, j, k])
                        piGrad =
                            kappaInv * MaInv2 * pEdgeR / rhou * (
                                (piR - piL) / dx +
                                met13EdgeR *
                                (piDDEdgeR - 4.0 * piDEdgeR + 3.0 * piEdgeR) *
                                0.5 / dz
                            )

                    else
                        piUEdgeR = 0.5 * (var.exner[i, j ,k+1] + var.exner[i+1, j, k+1])
                        piDEdgeR = 0.5 * (var.exner[i,j, k-1] + var.exner[i+1,j,k-1])
                        piGrad = kappaInv * MaInv2 * pEdgeR/rhou * ((piR - piL)/dx + met13EdgeR * (piUEdgeR - piDEdgeR)* 0.5/dz)
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
    end
end


function momentumPredictor_v!(semi, ode, mmp_mod, m)
    (; grid, cache, equations, met) = semi
    (; var, flux, rhoStrat, bvsStrat, pStrat, jac) = cache
    (; nx, ny, nz, lz, dx, dy, dz) = grid
    (; dt, alphaRK, betaRK, dMom, vOld, uOld, rhoOld) = ode
    (; kappaInv, MaInv2) = equations
    #periodic
    j0 = 0
    j1 = ny 
#    f_cor_nd[:] = 
    zBoundary = "solid_wall"
    if mmp_mod == "lhs"
    for k = 1:nz
        for j = j0:j1
            for i = 1:nx
                fluxDiff = compute_fluxDiff(flux.v, i, j, k, grid)

                jacEdgeF = 0.5f0 * (jac(i, j, k) + jac(i, j+1, k))
                fluxDiff = fluxDiff / jacEdgeF

                volForce = 0.0
                ## LHS
              #  vOld[i, j, k] = var.v[i, j, k]
              #  uC = 0.5 * (uOld[i, j, k] + uOld[i-1, j, k])
              #  uF = 0.5 * (uOld[i, j+1, k] + uOld[i-1, j+1, k])
                # TODO: Add Coriolis
              #  volForce =
              #     volForce -
              #     0.5 * (
              #         f_cor_nd[j] * (rhoOld[i, j, k] + rhoStrat[i, j, k]) * uC +
              #         f_or_nd[j+1] * (rhoOld[i, j+1, k] + rhoStrat[i, j+1, k]) * uF
              #     )
                
                F = -fluxDiff + volForce
                ## end LHS
                rhoM_1 = 0.5 * (rhoOld[i, j, k] + rhoOld[i, j+1, k])
                rhoM = 0.5 * (var.rho[i, j, k] + var.rho[i, j+1, k])

                rhoStratEdgeF = 0.5 * (rhoStrat[i, j, k] + rhoStrat[i, j+1, k])

                rhoM_1 += rhoStratEdgeF
                rhoM += rhoStratEdgeF

                vM_1 = var.v[i, j, k]
                momM_1 = rhoM_1 * vM_1

                # q(m-1) -> q(m)
                dMom[i, j, k, 2] = dt * F + alphaRK[m] * dMom[i, j, k, 2]

                # rhoV(m-1) -> rhoV(m)
                momM = momM_1 + betaRK[m] * dMom[i, j, k, 2]

                # calc v(m,*)
                vAst = momM / rhoM

                # vAst -> var
                var.v[i, j, k] = vAst

            end
        end
    end
elseif mmp_mod == "rhs"
    for k = 1:nz
        for j = j0:j1
            for i = 1:nx
                # Compute rhov
                rhov = 0.5 * (var.rho[i, j, k] + var.rho[i, j+1, k])
                rhoStratEdgeF = 0.5 * (rhoStrat[i, j, k] + rhoStrat[i, j+1, k])
                rhov += rhoStratEdgeF

                # Compute pi values
                piF = var.exner[i, j+1, k]
                piB = var.exner[i, j, k]

                # Compute values at cell edges
                pEdgeF = 0.5 * (pStrat[i, j, k] + pStrat[i, j+1, k])
                met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j+1, k, 2, 3))

                # Compute pressure gradient component
                if k == 1 && zBoundary == "solid_wall"
                    piUUEdgeF = 0.5 * (var.exner[i, j, k+2] + var.exner[i, j+1, k+2])
                    piUEdgeF = 0.5 * (var.exner[i, j, k+1] + var.exner[i, j+1, k+1])
                    piEdgeF = 0.5 * (var.exner[i, j, k] + var.exner[i, j+1, k])
                    piGrad =
                        kappaInv * MaInv2 * pEdgeF / rhov * (
                            (piF - piB) / dy +
                            met23EdgeF *
                            (-piUUEdgeF + 4.0 * piUEdgeF - 3.0 * piEdgeF) *
                            0.5 / dz
                        )

                elseif k == nz && zBoundary == "solid_wall"
                    piDDEdgeF = 0.5 * (var.exner[i, j, k-2] + var.exner[i, j+1, k-2])
                    piDEdgeF = 0.5 * (var.exner[i, j, k-1] + var.exner[i, j+1, k-1])
                    piEdgeF = 0.5 * (var.exner[i, j, k] + var.exner[i, j+1, k])
                    piGrad =
                        kappaInv * MaInv2 * pEdgeF / rhov * (
                            (piF - piB) / dy +
                            met23EdgeF *
                            (piDDEdgeF - 4.0 * piDEdgeF + 3.0 * piEdgeF) *
                            0.5 / dz
                        )

                else
                    piUEdgeF = 0.5 * (var.exner[i, j, k+1] + var.exner[i, j+1, k+1])
                    piDEdgeF = 0.5 * (var.exner[i, j, k-1] + var.exner[i, j+1, k-1])
                    piGrad =
                        kappaInv * MaInv2 * pEdgeF / rhov *
                        ((piF - piB) / dy + met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                end

                volfcy = 0.0

                # Compute vstar
                uhorx =
                    0.25 * (
                        var.u[i-1, j, k] +
                        var.u[i-1, j+1, k] +
                        var.u[i, j, k] +
                        var.u[i, j+1, k]
                    )
                vhory = var.v[i, j, k]
                # TODO: Coriolis to be added
              #  f_cor_v = 0.5 * (f_cor_nd[j] + f_cor_nd[j+1])

                # Coriolis force integrated on LHS
                vAst = vhory + dt * (-piGrad + volfcy / rhov)

                #if spongeLayer && sponge_uv
                 #   vAst -= dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j+1, k]) * vhory
                #end

                var.v[i, j, k] = vAst
            end
        end
    end
end
end #end


function momentumPredictor_w!(semi, ode, mmp_mod, m)
    (; grid, cache, equations, met) = semi
    (; var, flux, rhoStrat, bvsStrat, pStrat, jac, fluxDiffU, fluxDiffV) = cache
    (; nx, ny, nz, lz, dx, dy, dz) = grid
    (; dt, alphaRK, betaRK, dMom, vOld, uOld, rhoOld, rhopOld) = ode
    (; kappaInv, MaInv2, g_ndim) = equations
    k0 = 1
    k1 = nz - 1
    if mmp_mod == "lhs"
    for k = k0:k1
        for j = 1:ny
            for i = 1:nx
                # Compute vertical momentum flux divergence
                fR = flux.w[i, j, k, 1]
                fL = flux.w[i-1, j, k, 1]
                gF = flux.w[i, j, k, 2]
                gB = flux.w[i, j-1, k, 2]
                hU = flux.w[i, j, k, 3]
                hD = flux.w[i, j, k-1, 3]
                fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

                # Adjust Cartesian vertical momentum flux divergence
                jacEdgeU =
                    2.0 * jac(i, j, k) * jac(i, j, k+1) / (jac(i, j, k) + jac(i, j, k+1))
                fluxDiff /= jacEdgeU

                # Compute zonal momentum flux divergences
                for ll = 0:1, mm = 0:1
                    fR = flux.u[i-ll, j, k+mm, 1]
                    fL = flux.u[i-1-ll, j, k+mm, 1]
                    gF = flux.u[i-ll, j, k+mm, 2]
                    gB = flux.u[i-ll, j-1, k+mm, 2]
                    hU = flux.u[i-ll, j, k+mm, 3]
                    hD = flux.u[i-ll, j, k-1+mm, 3]

                    fluxDiffU[ll, mm] = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

                    jacEdgeR = 0.5 * (jac(i-ll, j, k+mm) + jac(i+1-ll, j, k+mm))
                    fluxDiffU[ll, mm] /= jacEdgeR
                end

                # Compute meridional momentum flux divergences
                for ll = 0:1, mm = 0:1
                    fR = flux.v[i, j-ll, k+mm, 1]
                    fL = flux.v[i-1, j-ll, k+mm, 1]
                    gF = flux.v[i, j-ll, k+mm, 2]
                    gB = flux.v[i, j-1-ll, k+mm, 2]
                    hU = flux.v[i, j-ll, k+mm, 3]
                    hD = flux.v[i, j-ll, k-1+mm, 3]

                    fluxDiffV[ll, mm] = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

                    jacEdgeF = 0.5 * (jac(i, j-ll, k+mm) + jac(i, j+1-ll, k+mm))
                    fluxDiffV[ll, mm] /= jacEdgeF
                end

                # Compute transformed vertical momentum flux divergence
                fluxDiff = trafo(
                    i,
                    j,
                    k,
                    fluxDiffU[0, 0],
                    fluxDiffU[0, 1],
                    fluxDiffU[1, 0],
                    fluxDiffU[1, 1],
                    fluxDiffV[0, 0],
                    fluxDiffV[0, 1],
                    fluxDiffV[1, 0],
                    fluxDiffV[1, 1],
                    fluxDiff,
                    "tfc", semi
                )

                # Explicit integration of Coriolis force in TFC
                if mmp_mod == "lhs"
                   # vC = 0.5 * (vOldTFC[i, j, k] + vOldTFC[i, j-1, k])
                   # vU = 0.5 * (vOldTFC[i, j, k+1] + vOldTFC[i, j-1, k+1])
                   # uC = 0.5 * (uOldTFC[i, j, k] + uOldTFC[i-1, j, k])
                   # uU = 0.5 * (uOldTFC[i, j, k+1] + uOldTFC[i-1, j, k+1])

                    # volForce +=
                    #     f_cor_nd[j] * (
                    #         jac[i, j, k+1] *
                    #         met[i, j, k, 1, 3] *
                    #         (rhoOld[i, j, k] + rhoStratTFC[i, j, k]) *
                    #         vC +
                    #         jac[i, j, k] *
                    #         met[i, j, k+1, 1, 3] *
                    #         (rhoOld[i, j, k+1] + rhoStratTFC[i, j, k+1]) *
                    #         vU
                    #     ) / (jac[i, j, k] + jac[i, j, k+1]) -
                    #     f_cor_nd[j] * (
                    #         jac[i, j, k+1] *
                    #         met[i, j, k, 2, 3] *
                    #         (rhoOld[i, j, k] + rhoStratTFC[i, j, k]) *
                    #         uC +
                    #         jac[i, j, k] *
                    #         met[i, j, k+1, 2, 3] *
                    #         (rhoOld[i, j, k+1] + rhoStratTFC[i, j, k+1]) *
                    #         uU
                    #     ) / (jac[i, j, k] + jac[i, j, k+1])
                end
                volForce = 0.0
                # Compute F(phi) = RHS
                if mmp_mod == "lhs"
                    F = -fluxDiff + volForce
                else
                    error("ERROR: wrong mmp_mod")
                end

                rhoM_1 =
                    (jac(i, j, k+1) * rhoOld[i, j, k] + jac(i, j, k) * rhoOld[i, j, k+1]) /
                    (jac(i, j, k) + jac(i, j, k+1))

                rhoM =
                    (
                        jac(i, j, k+1) * var.rho[i, j, k] +
                        jac(i, j, k) * var.rho[i, j, k+1]
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                rhoStratEdgeU =
                    (
                        jac(i, j, k+1) * rhoStrat[i, j, k] +
                        jac(i, j, k) * rhoStrat[i, j, k+1]
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                rhoM_1 += rhoStratEdgeU
                rhoM += rhoStratEdgeU

                wM_1 = var.w[i, j, k]
                momM_1 = rhoM_1 * wM_1

                # q(m-1) -> q(m)
                dMom[i, j, k, 3] = dt * F + alphaRK[m] * dMom[i, j, k, 3]

                # rhoW(m-1) -> rhoW(m)
                momM = momM_1 + betaRK[m] * dMom[i, j, k, 3]

                # Compute w(m,*)
                wAst = momM / rhoM

                # Assign wAst to var
                var.w[i, j, k] = wAst

            end
        end
    end
elseif mmp_mod == "rhs"

    for k = k0:k1
        for j = 1:ny
            for i = 1:nx
                rho000 = var.rho[i, j, k]
                rho001 = var.rho[i, j, k+1]

                rhow =
                    (
                        jac(i, j, k+1) * var.rho[i, j, k] +
                        jac(i, j, k) * var.rho[i, j, k+1]
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                rhoStratEdgeU =
                    (
                        jac(i, j, k+1) * rhoStrat[i, j, k] +
                        jac(i, j, k) * rhoStrat[i, j, k+1]
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                rho000 += rhoStrat[i, j, k]
                rho001 += rhoStrat[i, j, k+1]
                rhow += rhoStratEdgeU

                piU = var.exner[i, j, k+1]
                piD = var.exner[i, j, k]

                # Compute values at cell edges
                pEdgeU =
                    (
                        jac(i, j, k+1) * pStrat[i, j, k] +
                        jac(i, j, k) * pStrat[i, j, k+1]
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                met13EdgeU =
                    (
                        jac(i, j, k+1) * met(i, j, k, 1, 3) +
                        jac(i, j, k) * met(i, j, k+1, 1, 3)
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                met23EdgeU =
                    (
                        jac(i, j, k+1) * met(i, j, k, 2, 3) +
                        jac(i, j, k) * met(i, j, k+1, 2, 3)
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                met33EdgeU =
                    (
                        jac(i, j, k+1) * met(i, j, k, 3, 3) +
                        jac(i, j, k ) * met(i, j, k+1, 3, 3)
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                piREdgeU =
                    (
                        jac(i+1, j, k+1) * var.exner[i+1, j, k] +
                        jac(i+1, j, k) * var.exner[i+1, j, k+1]
                    ) / (jac(i+1, j, k) + jac(i+1, j, k+1))

                piLEdgeU =
                    (
                        jac(i-1, j, k+1) * var.exner[i-1, j, k] +
                        jac(i-1, j, k) * var.exner[i-1, j, k+1]
                    ) / (jac(i-1, j, k) + jac(i-1, j, k+1))

                piFEdgeU =
                    (
                        jac(i, j+1, k+1) * var.exner[i, j+1, k] +
                        jac(i, j+1, k) * var.exner[i, j+1, k+1]
                    ) / (jac(i, j+1, k) + jac(i, j+1, k+1))

                piBEdgeU =
                    (
                        jac(i, j-1, k+1) * var.exner[i, j-1, k] +
                        jac(i, j-1, k) * var.exner[i, j-1, k+1]
                    ) / (jac(i, j-1, k) + jac(i, j-1, k+1))

                # Compute pressure gradient component
                piGrad =
                    kappaInv * MaInv2 * pEdgeU / rhow * (
                        met13EdgeU * (piREdgeU - piLEdgeU) * 0.5 / dx +
                        met23EdgeU * (piFEdgeU - piBEdgeU) * 0.5 / dy +
                        met33EdgeU * (var.exner[i, j, k+1] - var.exner[i, j, k]) / dz
                    )

                volfcz = 0.0

                # wstar
                wvert = var.w[i, j, k]

                buoy =
                    -g_ndim * (
                        jac(i, j, k+1) * rhopOld[i, j, k] / rho000 / jac(i, j, k) +
                        jac(i, j, k) * rhopOld[i, j, k+1] / rho001 / jac(i, j, k+1)
                    ) / (jac(i, j, k) + jac(i, j, k+1))

                wAst = wvert + dt * (buoy - piGrad + volfcz / rhow)

              #  if spongeLayer
              #      wAst -=
              #          dt * (
              #              jac[i, j, k+1] * kr_sp_w_tfc[i, j, k] +
              #              jac[i, j, k] * kr_sp_w_tfc[i, j, k+1]
              #          ) / (jac[i, j, k] + jac[i, j, k+1]) * wvert
              #  end

                var.w[i, j, k] = wAst

            end
        end
    end
end

end ##
