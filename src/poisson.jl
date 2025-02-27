function val_PsIn(var, dt, opt, facray)
    # Calculates the matrix values for the pressure solver
    # The solver solves for dt * dp, hence no dt in the matrix elements

    # facray multiplies the Rayleigh-damping terms so that they are only
    # handled in the implicit time stepping (sponge and immersed boundary)

    # opt = "expl" =>
    # Pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = "impl" =>
    # Pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations

    # Local variables
    real = Float64
    fcscal, fcscal_u, fcscal_d::real
    AR, AB, AF, AD, AU, AC, ACH, ACV::real
    ARU, ARD, ALU, ALD, AFU, AFD, ABU, ABD::real
    AUU, ADD, ARUU, ARDD, ALUU, ALDD, AFUU, AFDD, ABUU, ABDD::real
    fcscal_r,
    fcscal_l,
    fcscal_f,
    fcscal_b,
    fcscal_ru,
    fcscal_rd,
    fcscal_lu,
    fcscal_ld,
    fcscal_fu,
    fcscal_fd,
    fcscal_bu,
    fcscal_bd::real
    fcscal_uu,
    fcscal_dd,
    fcscal_ruu,
    fcscal_rdd,
    fcscal_luu,
    fcscal_ldd,
    fcscal_fuu,
    fcscal_fdd,
    fcscal_buu,
    fcscal_bdd::real
    jacInv::real
    pEdgeRDiv, pEdgeLDiv, pEdgeFDiv, pEdgeBDiv, pEdgeUDiv, pEdgeDDiv::real
    pEdgeRGra, pEdgeLGra, pEdgeFGra, pEdgeBGra, pEdgeUGra, pEdgeDGra::real
    rhoEdgeR, rhoEdgeL, rhoEdgeF, rhoEdgeB, rhoEdgeU, rhoEdgeD::real
    met13EdgeR,
    met13EdgeL,
    met23EdgeF,
    met23EdgeB,
    met13EdgeU,
    met23EdgeU,
    met33EdgeU,
    met13EdgeD,
    met23EdgeD,
    met33EdgeD,
    met13UEdgeR,
    met13UEdgeL,
    met23UEdgeF,
    met23UEdgeB,
    met13DEdgeR,
    met13DEdgeL,
    met23DEdgeF,
    met23DEdgeB::real
    pUEdgeRGra,
    pUEdgeLGra,
    pUEdgeFGra,
    pUEdgeBGra,
    pDEdgeRGra,
    pDEdgeLGra,
    pDEdgeFGra,
    pDEdgeBGra::real
    rhoStratEdgeR,
    rhoStratEdgeL,
    rhoStratEdgeF,
    rhoStratEdgeB,
    rhoStratEdgeU,
    rhoStratEdgeD,
    rhoUEdgeR,
    rhoUEdgeL,
    rhoUEdgeF,
    rhoUEdgeB,
    rhoDEdgeR,
    rhoDEdgeL,
    rhoDEdgeF,
    rhoDEdgeB::real
    bvsStratEdgeU, bvsStratEdgeD::real
    facEdgeR,
    facEdgeL,
    facEdgeF,
    facEdgeB,
    facEdgeU,
    facEdgeD,
    facUEdgeR,
    facUEdgeL,
    facUEdgeF,
    facUEdgeB,
    facDEdgeR,
    facDEdgeL,
    facDEdgeF,
    facDEdgeB::real
    impHorEdgeR,
    impHorEdgeL,
    impHorEdgeF,
    impHorEdgeB,
    impHorUEdgeR,
    impHorUEdgeL,
    impHorUEdgeF,
    impHorUEdgeB,
    impHorDEdgeR,
    impHorDEdgeL,
    impHorDEdgeF,
    impHorDEdgeB,
    impVerEdgeU,
    impVerEdgeD::real
    gEdgeR,
    gEdgeL,
    gEdgeF,
    gEdgeB,
    gEdgeU,
    gEdgeD,
    gUEdgeR,
    gUEdgeL,
    gUEdgeF,
    gUEdgeB,
    gDEdgeR,
    gDEdgeL,
    gDEdgeF,
    gDEdgeB::real

    #---------------------------------
    # Loop over field
    #---------------------------------

    if opt == "expl"
        # Compute tensor elements for TFC.
        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    # Compute scaling factors.

                    ## PART 1
                    fcscal = sqrt(pStratTFC[i, j, k]^2 / rhoStratTFC[i, j, k])
                    fcscal_r = sqrt(pStratTFC[i+1, j, k]^2 / rhoStratTFC[i+1, j, k])
                    fcscal_l = sqrt(pStratTFC[i-1, j, k]^2 / rhoStratTFC[i-1, j, k])
                    fcscal_f = sqrt(pStratTFC[i, j+1, k]^2 / rhoStratTFC[i, j+1, k])
                    fcscal_b = sqrt(pStratTFC[i, j-1, k]^2 / rhoStratTFC[i, j-1, k])
                    fcscal_u = sqrt(pStratTFC[i, j, k+1]^2 / rhoStratTFC[i, j, k+1])
                    fcscal_d = sqrt(pStratTFC[i, j, k-1]^2 / rhoStratTFC[i, j, k-1])
                    fcscal_ru = sqrt(pStratTFC[i+1, j, k+1]^2 / rhoStratTFC[i+1, j, k+1])
                    fcscal_rd = sqrt(pStratTFC[i+1, j, k-1]^2 / rhoStratTFC[i+1, j, k-1])
                    fcscal_lu = sqrt(pStratTFC[i-1, j, k+1]^2 / rhoStratTFC[i-1, j, k+1])
                    fcscal_ld = sqrt(pStratTFC[i-1, j, k-1]^2 / rhoStratTFC[i-1, j, k-1])
                    fcscal_fu = sqrt(pStratTFC[i, j+1, k+1]^2 / rhoStratTFC[i, j+1, k+1])
                    fcscal_fd = sqrt(pStratTFC[i, j+1, k-1]^2 / rhoStratTFC[i, j+1, k-1])
                    fcscal_bu = sqrt(pStratTFC[i, j-1, k+1]^2 / rhoStratTFC[i, j-1, k+1])
                    fcscal_bd = sqrt(pStratTFC[i, j-1, k-1]^2 / rhoStratTFC[i, j-1, k-1])
                    fcscal_uu = sqrt(pStratTFC[i, j, k+2]^2 / rhoStratTFC[i, j, k+2])
                    fcscal_dd = sqrt(pStratTFC[i, j, k-2]^2 / rhoStratTFC[i, j, k-2])
                    fcscal_ruu = sqrt(pStratTFC[i+1, j, k+2]^2 / rhoStratTFC[i+1, j, k+2])
                    fcscal_rdd = sqrt(pStratTFC[i+1, j, k-2]^2 / rhoStratTFC[i+1, j, k-2])
                    fcscal_luu = sqrt(pStratTFC[i-1, j, k+2]^2 / rhoStratTFC[i-1, j, k+2])
                    fcscal_ldd = sqrt(pStratTFC[i-1, j, k-2]^2 / rhoStratTFC[i-1, j, k-2])
                    fcscal_fuu = sqrt(pStratTFC[i, j+1, k+2]^2 / rhoStratTFC[i, j+1, k+2])
                    fcscal_fdd = sqrt(pStratTFC[i, j+1, k-2]^2 / rhoStratTFC[i, j+1, k-2])
                    fcscal_buu = sqrt(pStratTFC[i, j-1, k+2]^2 / rhoStratTFC[i, j-1, k+2])
                    fcscal_bdd = sqrt(pStratTFC[i, j-1, k-2]^2 / rhoStratTFC[i, j-1, k-2])

                    ## PART 2

                    # Compute inverse Jacobian
                    jacInv = 1.0 / jac(i, j, k)

                    # Compute P coefficients (divergence)
                    pEdgeRDiv =
                        0.5 * (
                            jac(i, j, k) * pStratTFC[i, j, k] +
                            jac(i + 1, j, k) * pStratTFC[i+1, j, k]
                        )
                    pEdgeLDiv =
                        0.5 * (
                            jac(i, j, k) * pStratTFC[i, j, k] +
                            jac(i - 1, j, k) * pStratTFC[i-1, j, k]
                        )
                    pEdgeFDiv =
                        0.5 * (
                            jac(i, j, k) * pStratTFC[i, j, k] +
                            jac(i, j + 1, k) * pStratTFC[i, j+1, k]
                        )
                    pEdgeBDiv =
                        0.5 * (
                            jac(i, j, k) * pStratTFC[i, j, k] +
                            jac(i, j - 1, k) * pStratTFC[i, j-1, k]
                        )
                    pEdgeUDiv =
                        jac(i, j, k) *
                        jac(i, j, k + 1) *
                        (pStratTFC[i, j, k] + pStratTFC[i, j, k+1]) /
                        (jac(i, j, k) + jac(i, j, k + 1))
                    pEdgeDDiv =
                        jac(i, j, k) *
                        jac(i, j, k - 1) *
                        (pStratTFC[i, j, k] + pStratTFC[i, j, k-1]) /
                        (jac(i, j, k) + jac(i, j, k - 1))

                    # Compute P coefficients (pressure gradient)
                    pEdgeRGra = 0.5 * (pStratTFC[i, j, k] + pStratTFC[i+1, j, k])
                    pEdgeLGra = 0.5 * (pStratTFC[i, j, k] + pStratTFC[i-1, j, k])
                    pEdgeFGra = 0.5 * (pStratTFC[i, j, k] + pStratTFC[i, j+1, k])
                    pEdgeBGra = 0.5 * (pStratTFC[i, j, k] + pStratTFC[i, j-1, k])
                    pEdgeUGra =
                        (
                            jac(i, j, k + 1) * pStratTFC[i, j, k] +
                            jac(i, j, k) * pStratTFC[i, j, k+1]
                        ) / (jac(i, j, k) + jac(i, j, k + 1))
                    pEdgeDGra =
                        (
                            jac(i, j, k - 1) * pStratTFC[i, j, k] +
                            jac(i, j, k) * pStratTFC[i, j, k-1]
                        ) / (jac(i, j, k) + jac(i, j, k - 1))

                    # Compute density coefficients
                    rhoEdgeR =
                        0.5 * (
                            var.rho[i, j, k] +
                            var.rho[i+1, j, k] +
                            rhoStratTFC[i, j, k] +
                            rhoStratTFC[i+1, j, k]
                        )
                    rhoEdgeL =
                        0.5 * (
                            var.rho[i, j, k] +
                            var.rho[i-1, j, k] +
                            rhoStratTFC[i, j, k] +
                            rhoStratTFC[i-1, j, k]
                        )
                    rhoEdgeF =
                        0.5 * (
                            var.rho[i, j, k] +
                            var.rho[i, j+1, k] +
                            rhoStratTFC[i, j, k] +
                            rhoStratTFC[i, j+1, k]
                        )
                    rhoEdgeB =
                        0.5 * (
                            var.rho[i, j, k] +
                            var.rho[i, j-1, k] +
                            rhoStratTFC[i, j, k] +
                            rhoStratTFC[i, j-1, k]
                        )
                    rhoEdgeU =
                        (
                            jac(i, j, k + 1) * (var.rho[i, j, k] + rhoStratTFC[i, j, k]) +
                            jac(i, j, k) * (var.rho[i, j, k+1] + rhoStratTFC[i, j, k+1])
                        ) / (jac(i, j, k) + jac(i, j, k + 1))
                    rhoEdgeD =
                        (
                            jac(i, j, k - 1) * (var.rho[i, j, k] + rhoStratTFC[i, j, k]) +
                            jac(i, j, k) * (var.rho[i, j, k-1] + rhoStratTFC[i, j, k-1])
                        ) / (jac(i, j, k) + jac(i, j, k - 1))

                    # Interpolate metric-tensor elements
                    met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                    met13EdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
                    met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                    met23EdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
                    met13EdgeU =
                        (
                            jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                            jac(i, j, k) * met(i, j, k + 1, 1, 3)
                        ) / (jac(i, j, k) + jac(i, j, k + 1))
                    met23EdgeU =
                        (
                            jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                            jac(i, j, k) * met(i, j, k + 1, 2, 3)
                        ) / (jac(i, j, k) + jac(i, j, k + 1))
                    met33EdgeU =
                        (
                            jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                            jac(i, j, k) * met(i, j, k + 1, 3, 3)
                        ) / (jac(i, j, k) + jac(i, j, k + 1))
                    met13EdgeD =
                        (
                            jac(i, j, k - 1) * met(i, j, k, 1, 3) +
                            jac(i, j, k) * met(i, j, k - 1, 1, 3)
                        ) / (jac(i, j, k) + jac(i, j, k - 1))
                    met23EdgeD =
                        (
                            jac(i, j, k - 1) * met(i, j, k, 2, 3) +
                            jac(i, j, k) * met(i, j, k - 1, 2, 3)
                        ) / (jac(i, j, k) + jac(i, j, k - 1))
                    met33EdgeD =
                        (
                            jac(i, j, k - 1) * met(i, j, k, 3, 3) +
                            jac(i, j, k) * met(i, j, k - 1, 3, 3)
                        ) / (jac(i, j, k) + jac(i, j, k - 1))

                    # --------------------- A(i,j,k) ---------------------

                    if k == 1 && zBoundary == "solid_wall"
                        AC =
                            -jacInv / dx * (
                                pEdgeRDiv / rhoEdgeR *
                                pEdgeRGra *
                                (1.0 / dx + 0.75 * met13EdgeR / dz) +
                                pEdgeLDiv / rhoEdgeL *
                                pEdgeLGra *
                                (1.0 / dx - 0.75 * met13EdgeL / dz)
                            ) -
                            jacInv / dy * (
                                pEdgeFDiv / rhoEdgeF *
                                pEdgeFGra *
                                (1.0 / dy + 0.75 * met23EdgeF / dz) +
                                pEdgeBDiv / rhoEdgeB *
                                pEdgeBGra *
                                (1.0 / dy - 0.75 * met23EdgeB / dz)
                            ) -
                            jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU / dz
                    elseif k == nz && zBoundary == "solid_wall"
                        AC =
                            -jacInv / dx * (
                                pEdgeRDiv / rhoEdgeR *
                                pEdgeRGra *
                                (1.0 / dx - 0.75 * met13EdgeR / dz) +
                                pEdgeLDiv / rhoEdgeL *
                                pEdgeLGra *
                                (1.0 / dx + 0.75 * met13EdgeL / dz)
                            ) -
                            jacInv / dy * (
                                pEdgeFDiv / rhoEdgeF *
                                pEdgeFGra *
                                (1.0 / dy - 0.75 * met23EdgeF / dz) +
                                pEdgeBDiv / rhoEdgeB *
                                pEdgeBGra *
                                (1.0 / dy + 0.75 * met23EdgeB / dz)
                            ) -
                            jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD / dz
                    else
                        AC =
                            -jacInv / dx * (
                                pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx +
                                pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx
                            ) -
                            jacInv / dy * (
                                pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy +
                                pEdgeBDiv / rhoEdgeB * pEdgeBGra / dy
                            ) -
                            jacInv / dz * (
                                pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU / dz +
                                pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD / dz
                            )
                    end

                    ## Part 3



                end
            end
        end
    end
end
