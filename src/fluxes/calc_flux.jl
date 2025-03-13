
function compute_fluxes!(model)
    @trixi_timeit timer() "Compute fluxes" begin
    #! format: noindent
    reconstruction!(model)
    rhoFlux!(model)
    rhopFlux!(model)
    momemtumFlux!(model)
    end # timer
end

# TODO: reconstruct_rho! could probably be called in rhoFlux!

function rhoFlux!(model)
    (; nx, ny, nz) = model.domain

    rhoFlux = model.fluxes.rho # mass flux
    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    rhoStrat = model.atmosphere.rhostrattfc
    rhoTilde = model.variables.reconstructed.rho
    jac = model.grid.jac

    # TODO: these could be moved to new functions
    # zonal rho flux
    for kz in 1:nz
        for jy in 1:ny
            for ix in 0:nx
                rhoStratEdgeR = 0.5 * (rhoStrat[ix, jy, kz] + rhoStrat[ix+1, jy, kz])
                pEdgeR = 0.5 * (pStrat[ix, jy, kz] + pStrat[ix+1, jy, kz])
                rhoR = rhoTilde[ix+1, jy, kz, 1, 0] + rhoStratEdgeR / pEdgeR
                rhoL = rhoTilde[ix, jy, kz, 1, 1] + rhoStratEdgeR / pEdgeR

                pEdgeR = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                                jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz])
                uSurf = pEdgeR * u[ix, jy, kz]
                # TODO: make sure this is the correct u (see vara in fluxes.f90)

                fRho = flux_muscl(uSurf, rhoL, rhoR)

                rhoFlux[ix, jy, kz, 1] = fRho
            end
        end
    end

    # meridional rho flux
    for kz in 1:nz
        for jy in 0:ny
            for ix in 1:nx
                rhoStratEdgeF = 0.5 * (rhoStrat[ix, jy, kz] + rhoStrat[ix, jy+1, kz])
                pEdgeF = 0.5 * (pStrat[ix, jy, kz] + pStrat[ix, jy+1, kz])
                rhoF = rhoTilde[ix, jy+1, kz, 2, 0] + rhoStratEdgeF / pEdgeF
                rhoB = rhoTilde[ix, jy, kz, 2, 1] + rhoStratEdgeF / pEdgeF

                pEdgeF = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                                jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz])
                vSurf = pEdgeF * v[ix, jy, kz]

                gRho = flux_muscl(vSurf, rhoB, rhoF)

                rhoFlux[ix, jy, kz, 2] = gRho
            end
        end
    end

    # vertical rho flux
    for kz in 0:nz
        for jy in 1:ny
            for ix in 1:nx
                rhoStratEdgeU = (jac(ix, jy, kz + 1) * rhoStrat[ix, jy, kz] +
                                 jac(ix, jy, kz) * rhoStrat[ix, jy, kz+1]) /
                                (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pEdgeU = (jac(ix, jy, kz + 1) * pStrat[ix, jy, kz] +
                          jac(ix, jy, kz) * pStrat[ix, jy, kz+1]) /
                         (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                rhoU = rhoTilde[ix, jy, kz+1, 3, 0] + rhoStratEdgeU / pEdgeU
                rhoD = rhoTilde[ix, jy, kz, 3, 1] + rhoStratEdgeU / pEdgeU

                pEdgeU = jac(ix, jy, kz) *
                         jac(ix, jy, kz + 1) *
                         (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                         (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                wSurf = pEdgeU * w[ix, jy, kz]

                hRho = flux_muscl(wSurf, rhoD, rhoU)

                rhoFlux[ix, jy, kz, 3] = hRho
            end
        end
    end
end

function rhopFlux!(model)
    (; nx, ny, nz) = model.domain

    rhopFlux = model.fluxes.rhop # mass flux
    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    rhopTilde = model.variables.reconstructed.rhop
    jac = model.grid.jac

    for kz in 1:nz
        for jy in 1:ny
            for ix in 0:nx
                rhopR = rhopTilde[ix+1, jy, kz, 1, 0]
                rhopL = rhopTilde[ix, jy, kz, 1, 1]

                pEdgeR = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                                jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz])
                uSurf = pEdgeR * u[ix, jy, kz]

                fRhop = flux_muscl(uSurf, rhopL, rhopR)

                rhopFlux[ix, jy, kz, 1] = fRhop
            end
        end
    end

    # meridional rho flux
    for kz in 1:nz
        for jy in 0:ny
            for ix in 1:nx
                rhopF = rhopTilde[ix, jy+1, kz, 2, 0]
                rhopB = rhopTilde[ix, jy, kz, 2, 1]

                pEdgeF = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                                jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz])
                vSurf = pEdgeF * v[ix, jy, kz]

                gRhop = flux_muscl(vSurf, rhopB, rhopF)

                rhopFlux[ix, jy, kz, 2] = gRhop
            end
        end
    end

    # vertical rho flux
    for kz in 0:nz
        for jy in 1:ny
            for ix in 1:nx
                rhopU = rhopTilde[ix, jy, kz+1, 3, 0]
                rhopD = rhopTilde[ix, jy, kz, 3, 1]

                pEdgeU = jac(ix, jy, kz) *
                         jac(ix, jy, kz + 1) *
                         (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                         (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                wSurf = pEdgeU * w[ix, jy, kz]

                hRhop = flux_muscl(wSurf, rhopD, rhopU)

                rhopFlux[ix, jy, kz, 3] = hRhop
            end
        end
    end
end

function momemtumFlux!(semi)
    rhouFlux!(semi) # (rho*u)*u, (rho*u)*v, (rho*u)*w
    rhovFlux!(semi) # (rho*v)*u, (rho*v)*v, (rho*v)*w
    rhowFlux!(semi) # (rho*w)*u, (rho*w)*v, (rho*w)*w
end

function rhouFlux!(semi)
    rhouuFlux!(semi) # calculate (rho*u)*u
    rhouvFlux!(semi) # calculate (rho*u)*v
    rhouwFlux!(semi) # calculate (rho*u)*w
end

function rhovFlux!(semi)
    rhovuFlux!(semi) # calculate (rho*v)*u
    rhovvFlux!(semi) # calculate (rho*v)*v
    rhovwFlux!(semi) # calculate (rho*v)*w
end

function rhowFlux!(semi)
    rhowuFlux!(semi) # calculate (rho*w)*u
    rhowvFlux!(semi) # calculate (rho*w)*v
    rhowwFlux!(semi) # calculate (rho*w)*w
end

function rhouuFlux!(model) # (rho*u)*u
    uFlux = model.fluxes.u # mass flux
    (; nx, ny, nz) = model.domain

    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    uTilde = model.variables.reconstructed.u
    jac = model.grid.jac

    # Flux fRhoU
    for kz in 1:nz
        for jy in 1:ny
            for ix in -1:nx
                uR = uTilde[ix+1, jy, kz, 1, 0]
                uL = uTilde[ix, jy, kz, 1, 1]

                pEdge = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                               jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz])
                pREdge = 0.5 * (jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz] +
                                jac(ix + 2, jy, kz) * pStrat[ix+2, jy, kz])
                uSurf = 0.5 * (pEdge * u[ix, jy, kz] + pREdge * u[ix+1, jy, kz])

                fRhoU = flux_muscl(uSurf, uL, uR)

                uFlux[ix, jy, kz, 1] = fRhoU
            end
        end
    end
end

function rhouvFlux!(model) # (rho*u)*v

    uFlux = model.fluxes.u # mass flux
    (; nx, ny, nz) = model.domain

    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    uTilde = model.variables.reconstructed.u
    jac = model.grid.jac
    # Flux gRhoU
    for kz in 1:nz
        for jy in 0:ny
            for ix in 0:nx
                uF = uTilde[ix, jy+1, kz, 2, 0]
                uB = uTilde[ix, jy, kz, 2, 1]

                pEdge = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                               jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz])
                pREdge = 0.5 * (jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz] +
                                jac(ix + 1, jy + 1, kz) * pStrat[ix+1, jy+1, kz])
                vSurf = 0.5 * (pEdge * v[ix, jy, kz] + pREdge * v[ix+1, jy, kz])

                gRhoU = flux_muscl(vSurf, uB, uF)

                uFlux[ix, jy, kz, 2] = gRhoU
            end
        end
    end
end

function rhouwFlux!(model) # (rho*u)*w

    uFlux = model.fluxes.u # mass flux
    (; nx, ny, nz) = model.domain

    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    uTilde = model.variables.reconstructed.u
    jac = model.grid.jac

    # Flux hRhoU
    for kz in 0:nz
        for jy in 1:ny
            for ix in 0:nx
                uU = uTilde[ix, jy, kz+1, 3, 0]
                uD = uTilde[ix, jy, kz, 3, 1]

                pEdge = jac(ix, jy, kz) *
                        jac(ix, jy, kz + 1) *
                        (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                        (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pREdge = jac(ix + 1, jy, kz) *
                         jac(ix + 1, jy, kz + 1) *
                         (pStrat[ix+1, jy, kz] + pStrat[ix+1, jy, kz+1]) /
                         (jac(ix + 1, jy, kz) + jac(ix + 1, jy, kz + 1))
                wSurf = 0.5 * (pEdge * w[ix, jy, kz] + pREdge * w[ix+1, jy, kz])

                hRhoU = flux_muscl(wSurf, uD, uU)

                uFlux[ix, jy, kz, 3] = hRhoU
            end
        end
    end
end

function rhovuFlux!(model)
    (; nx, ny, nz) = model.domain

    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    vTilde = model.variables.reconstructed.v
    vFlux = model.fluxes.v # mass flux
    jac = model.grid.jac

    # Flux fRhoV
    for kz in 1:nz
        for jy in 0:ny
            for ix in 0:nx
                vR = vTilde[ix+1, jy, kz, 1, 0]
                vL = vTilde[ix, jy, kz, 1, 1]

                pEdge = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                               jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz])
                pREdge = 0.5 * (jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz] +
                                jac(ix + 1, jy + 1, kz) * pStrat[ix+1, jy+1, kz])
                uSurf = 0.5 * (pEdge * u[ix, jy, kz] + pREdge * u[ix, jy+1, kz])

                fRhoV = flux_muscl(uSurf, vL, vR)

                vFlux[ix, jy, kz, 1] = fRhoV
            end
        end
    end
end

function rhovvFlux!(model) # (rho*v)*v
    (; nx, ny, nz) = model.domain

    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    vTilde = model.variables.reconstructed.v
    vFlux = model.fluxes.v # mass flux
    jac = model.grid.jac
    # Flux gRhoV
    for kz in 1:nz
        for jy in -1:ny
            for ix in 1:nx
                vF = vTilde[ix, jy+1, kz, 2, 0]
                vB = vTilde[ix, jy, kz, 2, 1]

                pEdge = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                               jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz])
                pREdge = 0.5 * (jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz] +
                                jac(ix, jy + 2, kz) * pStrat[ix, jy+2, kz])
                vSurf = 0.5 * (pEdge * v[ix, jy, kz] + pREdge * v[ix, jy+1, kz])

                gRhoV = flux_muscl(vSurf, vB, vF)

                vFlux[ix, jy, kz, 2] = gRhoV
            end
        end
    end
end

function rhovwFlux!(model)
    (; nx, ny, nz) = model.domain
    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    vTilde = model.variables.reconstructed.v
    vFlux = model.fluxes.v # mass flux
    jac = model.grid.jac
    # Flux hRhoV
    for kz in 0:nz
        for jy in 0:ny
            for ix in 1:nx
                vU = vTilde[ix, jy, kz+1, 3, 0]
                vD = vTilde[ix, jy, kz, 3, 1]

                pEdge = jac(ix, jy, kz) *
                        jac(ix, jy, kz + 1) *
                        (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                        (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pREdge = jac(ix, jy + 1, kz) *
                         jac(ix, jy + 1, kz + 1) *
                         (pStrat[ix, jy+1, kz] + pStrat[ix, jy+1, kz+1]) /
                         (jac(ix, jy + 1, kz) + jac(ix, jy + 1, kz + 1))
                wSurf = 0.5 * (pEdge * w[ix, jy, kz] + pREdge * w[ix, jy+1, kz])
                hRhoV = flux_muscl(wSurf, vD, vU)

                vFlux[ix, jy, kz, 3] = hRhoV
            end
        end
    end
end

function rhowuFlux!(model) # (rho*w)*u
    (; nx, ny, nz) = model.domain
    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    wTilde = model.variables.reconstructed.w
    wFlux = model.fluxes.w # mass flux
    jac = model.grid.jac

    # Flux fRhoW
    for kz in 0:nz
        for jy in 1:ny
            for ix in 0:nx
                wR = wTilde[ix+1, jy, kz, 1, 0]
                wL = wTilde[ix, jy, kz, 1, 1]

                pEdge = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                               jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz])
                pREdge = 0.5 * (jac(ix, jy, kz + 1) * pStrat[ix, jy, kz+1] +
                                jac(ix + 1, jy, kz + 1) * pStrat[ix+1, jy, kz+1])
                uSurf = ((jac[ix, jy, kz+1] + jac[ix+1, jy, kz+1]) * pEdge *
                         u[ix, jy, kz] +
                         (jac[ix, jy, kz] + jac[ix+1, jy, kz]) * pREdge *
                         u[ix, jy, kz+1]) / (jac[ix, jy, kz] +
                                             jac[ix+1, jy, kz] +
                                             jac[ix, jy, kz+1] +
                                             jac[ix+1, jy, kz+1])

                fRhoW = flux_muscl(uSurf, wL, wR)

                wFlux[ix, jy, kz, 1] = fRhoW
            end
        end
    end
end

function rhowvFlux!(model)
    (; nx, ny, nz) = model.domain
    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    wTilde = model.variables.reconstructed.w
    wFlux = model.fluxes.w # mass flux
    jac = model.grid.jac

    # Flux gRhoW
    for kz in 0:nz
        for jy in 0:ny
            for ix in 1:nx
                wF = wTilde[ix, jy+1, kz, 2, 0]
                wB = wTilde[ix, jy, kz, 2, 1]

                pEdge = 0.5 * (jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                               jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz])
                pREdge = 0.5 * (jac(ix, jy, kz + 1) * pStrat[ix, jy, kz+1] +
                                jac(ix, jy + 1, kz + 1) * pStrat[ix, jy+1, kz+1])
                vSurf = ((jac(ix, jy, kz + 1) + jac(ix, jy + 1, kz + 1)) *
                         pEdge *
                         v[ix, jy, kz] +
                         (jac(ix, jy, kz) + jac(ix, jy + 1, kz)) * pREdge *
                         v[ix, jy, kz+1]) / (jac(ix, jy, kz) +
                                             jac(ix, jy + 1, kz) +
                                             jac(ix, jy, kz + 1) +
                                             jac(ix, jy + 1, kz + 1))

                gRhoW = flux_muscl(vSurf, wB, wF)

                wFlux[ix, jy, kz, 2] = gRhoW
            end
        end
    end
end

function rhowwFlux!(model)
    (; nx, ny, nz) = model.domain
    (; u, v, w) = model.variables.prognostic_fields_0
    pStrat = model.atmosphere.pstrattfc
    wTilde = model.variables.reconstructed.w
    wFlux = model.fluxes.w # mass flux
    jac = model.grid.jac

    # Flux hRhoW
    for kz in -1:nz
        for jy in 1:ny
            for ix in 1:nx
                wU = wTilde[ix, jy, kz+1, 3, 0]
                wD = wTilde[ix, jy, kz, 3, 1]

                pEdge = jac(ix, jy, kz) *
                        jac(ix, jy, kz + 1) *
                        (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                        (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pREdge = jac(ix, jy, kz + 1) *
                         jac(ix, jy, kz + 2) *
                         (pStrat[ix, jy, kz+1] + pStrat[ix, jy, kz+2]) /
                         (jac(ix, jy, kz + 1) + jac(ix, jy, kz + 2))
                wSurf = 0.5 * (pEdge * w[ix, jy, kz] + pREdge * w[ix, jy, kz+1])

                hRhoW = flux_muscl(wSurf, wD, wU)

                wFlux[ix, jy, kz, 3] = hRhoW
            end
        end
    end
end

function flux_muscl(wind, phiUp, phiDown)
    if (wind > 0.0)
        return wind * phiUp
    else
        return wind * phiDown
    end
end
