#--------- RECONSTRUCTION ---------#

# TODO: check start and end of loops
function reconstruction!(semi)

    reconstruct_rho!(semi)
    reconstruct_rhop!(semi)
    reconstruct_u!(semi)
    reconstruct_v!(semi)
    reconstruct_w!(semi)

end


function reconstruct_rho!(semi)

    (; cache, grid) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, rhoTilde, pStrat) = cache

    rho = var.rho
    sizeX, sizeY, sizeZ = size(rho)

    rhoBar_ = zeros(sizeX, sizeY, sizeZ)
    rhoBar = OffsetArray(rhoBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    for ix = -nbx:nx+nbx
        for jy = -nby:ny+nby
            for kz = 0:nz+1
                # TODO: error message for pStrat
                rhoBar[ix, jy, kz] = rho[ix, jy, kz] / pStrat[ix, jy, kz]
            end
        end
    end
    muscl_reconstruct3D!(rhoBar, nxx, nyy, nzz, rhoTilde)
end

function reconstruct_rhop!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, rhopTilde, pStrat) = cache

    rhop = var.rhop
    sizeX, sizeY, sizeZ = size(rhop)

    rhopBar_ = zeros(sizeX, sizeY, sizeZ)
    rhopBar = OffsetArray(rhopBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    for ix = -nbx:nx+nbx
        for jy = -nby:ny+nby
            for kz = 0:nz+1
                rhopBar[ix, jy, kz] = rhop[ix, jy, kz] / pStrat[ix, jy, kz]
            end
        end
    end
    muscl_reconstruct3D!(rhopBar, nxx, nyy, nzz, rhopTilde)
    # TODO: keep limiterType1 ?
end

function reconstruct_u!(semi)
    # Compute rho*u/P for reconstruction
    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, uTilde, rhoStrat, pStrat) = cache

    rho = var.rho
    u = var.u
    sizeX, sizeY, sizeZ = size(u)

    uBar_ = zeros(sizeX, sizeY, sizeZ)
    uBar = OffsetArray(uBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    # rho*u/P for reconstruction
    for ix = -nbx:nx+nbx-1
        for jy = -nby:ny+nby
            for kz = 0:nz+1
                rhoEdge =
                    0.5 * (
                        rho[ix, jy, kz] +
                        rho[ix+1, jy, kz] +
                        rhoStrat[ix, jy, kz] +
                        rhoStrat[ix+1, jy, kz]
                    )
                pEdge = 0.5 * (pStrat[ix, jy, kz] + pStrat[ix+1, jy, kz])
                uBar[ix, jy, kz] = u[ix, jy, kz] * rhoEdge / pEdge
            end
        end
    end

    muscl_reconstruct3D!(uBar, nxx, nyy, nzz, uTilde)
end

function reconstruct_v!(semi)

    # Compute rho*v/P for reconstruction
    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, vTilde, rhoStrat, pStrat) = cache

    rho = var.rho
    v = var.v
    sizeX, sizeY, sizeZ = size(v)

    vBar_ = zeros(sizeX, sizeY, sizeZ)
    vBar = OffsetArray(vBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    # rho*v/P for reconstruction
    for ix = -nbx:nx+nbx
        for jy = -nby:ny+nby-1
            for kz = 0:nz+1
                rhoEdge =
                    0.5 * (
                        rho[ix, jy, kz] +
                        rhoStrat[ix, jy, kz] +
                        rho[ix, jy+1, kz] +
                        rhoStrat[ix, jy+1, kz]
                    )
                pEdge = 0.5 * (pStrat[ix, jy, kz] + pStrat[ix, jy+1, kz])
                vBar[ix, jy, kz] = v[ix, jy, kz] * rhoEdge / pEdge
            end
        end
    end

    muscl_reconstruct3D!(vBar, nxx, nyy, nzz, vTilde)
end

function reconstruct_w!(semi)

    # Compute rho*v/P for reconstruction
    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, lz) = grid
    (; var, wTilde, rhoStrat, pStrat, jac) = cache

    rho = var.rho
    w = var.w
    sizeX, sizeY, sizeZ = size(w)

    wBar_ = zeros(sizeX, sizeY, sizeZ)
    wBar = OffsetArray(wBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    wBar[:, :, 0:nz+1] .= w[:, :, 0:nz+1]
    # TODO: vertWind doesn't seem to work!!!
    for ix = 1:nx
        for jy = 1:ny
            for kz = 0:nz+1
                wBar[ix, jy, kz] = vertWind(ix, jy, kz, semi)
            end
        end
    end


    setHalosOfField!(wBar, semi)
    for ix = -nbx:nx+nbx
        for jy = -nby:ny+nby
            for kz = 0:nz+1
                rhoEdge =
                    (
                        jac(ix, jy, kz + 1) * (rho[ix, jy, kz] + rhoStrat[ix, jy, kz]) +
                        jac(ix, jy, kz) * (rho[ix, jy, kz+1] + rhoStrat[ix, jy, kz+1])
                    ) / (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pEdge =
                    (
                        jac(ix, jy, kz + 1) * pStrat[ix, jy, kz] +
                        jac(ix, jy, kz) * pStrat[ix, jy, kz+1]
                    ) / (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                wBar[ix, jy, kz] = wBar[ix, jy, kz] * rhoEdge / pEdge
            end
        end
    end
    muscl_reconstruct3D!(wBar, nxx, nyy, nzz, wTilde; debug = 1)

end

function muscl_reconstruct3D!(varBar, sizeX, sizeY, sizeZ, varTilde; debug = 0)

    # TODO: changed indexing. Check if output is correct!

    phiX = zeros(sizeX)
    phiY = zeros(sizeY)
    phiZ = zeros(sizeZ)

    phiTildeX = OffsetArray(zeros(sizeX, 2), 1:sizeX, 0:1)
    phiTildeY = OffsetArray(zeros(sizeY, 2), 1:sizeY, 0:1)
    phiTildeZ = OffsetArray(zeros(sizeZ, 2), 1:sizeZ, 0:1)

    for kz = 2:sizeZ-1
        for jy = 2:sizeY-1

            phiX[:] .= varBar.parent[:, jy, kz]
            muscle_reconstruct1D_mcvariant!(phiX, sizeX, phiTildeX)
            varTilde.parent[:, jy, kz, 1, :] .= phiTildeX.parent
        end
    end

    for kz = 2:sizeZ-1
        for ix = 2:sizeX-1
            phiY .= varBar.parent[ix, :, kz]
            muscle_reconstruct1D_mcvariant!(phiY, sizeY, phiTildeY)
            varTilde.parent[ix, :, kz, 2, :] .= phiTildeY.parent
        end
    end

    for jy = 2:sizeY-1
        for ix = 2:sizeX-1
            phiZ .= varBar.parent[ix, jy, :]
            muscle_reconstruct1D_mcvariant!(phiZ, sizeZ, phiTildeZ)
            varTilde.parent[ix, jy, :, 3, :] .= phiTildeZ.parent
        end
    end
end

function muscle_reconstruct1D_mcvariant!(phi, phiSize, phiTilde)

    for i = 2:phiSize-1
        deltaL = phi[i] - phi[i-1]
        deltaR = phi[i+1] - phi[i]

        if (deltaR == 0.0)
            phiTilde[i, 1] = phi[i] # 2: phiTilde_R
            phiTilde[i, 0] = phi[i] # 1: phiTilde_L
        else
            if (deltaL == 0.0)
                theta = deltaL / deltaR
                s = (2.0 + theta) / 3.0
                sigmaL = max(0.0, min(2 * theta, s, 2.0))

                phiTilde[i, 1] = phi[i] + 0.5 * sigmaL * deltaR
                phiTilde[i, 0] = phi[i]
            else
                theta = deltaL / deltaR
                s = (2.0 + theta) / 3.0
                sigmaL = max(0.0, min(2 * theta, s, 2.0))
                s = (2.0 + 1.0 / theta) / 3.0
                sigmaR = max(0.0, min(2 / theta, s, 2.0))

                phiTilde[i, 1] = phi[i] + 0.5 * sigmaL * deltaR
                phiTilde[i, 0] = phi[i] - 0.5 * sigmaR * deltaL
            end
        end

    end



end

# ------------------------ FLUX CALCULATION ------------------------
function compute_fluxes!(semi)

    reconstruction!(semi)
    rhoFlux!(semi)
    rhopFlux!(semi)
    momemtumFlux!(semi)

end

# TODO: reconstruct_rho! could probably be called in rhoFlux!

function rhoFlux!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, lz) = grid
    (; var, var0, rhoTilde, rhoStrat, pStrat, flux, jac) = cache

    rhoFlux = flux.rho # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # TODO: these could be moved to new functions
    # zonal rho flux
    for kz = 1:nz
        for jy = 1:ny
            for ix = 0:nx
                rhoStratEdgeR = 0.5 * (rhoStrat[ix, jy, kz] + rhoStrat[ix+1, jy, kz])
                pEdgeR = 0.5 * (pStrat[ix, jy, kz] + pStrat[ix+1, jy, kz])
                rhoR = rhoTilde[ix+1, jy, kz, 1, 0] + rhoStratEdgeR / pEdgeR
                rhoL = rhoTilde[ix, jy, kz, 1, 1] + rhoStratEdgeR / pEdgeR

                pEdgeR =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz]
                    )
                uSurf = pEdgeR * u[ix, jy, kz]
                # TODO: make sure this is the correct u (see vara in fluxes.f90)

                fRho = flux_muscl(uSurf, rhoL, rhoR)

                rhoFlux[ix, jy, kz, 1] = fRho
            end
        end
    end

    # meridional rho flux
    for kz = 1:nz
        for jy = 0:ny
            for ix = 1:nx
                rhoStratEdgeF = 0.5 * (rhoStrat[ix, jy, kz] + rhoStrat[ix, jy+1, kz])
                pEdgeF = 0.5 * (pStrat[ix, jy, kz] + pStrat[ix, jy+1, kz])
                rhoF = rhoTilde[ix, jy+1, kz, 2, 0] + rhoStratEdgeF / pEdgeF
                rhoB = rhoTilde[ix, jy, kz, 2, 1] + rhoStratEdgeF / pEdgeF

                pEdgeF =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz]
                    )
                vSurf = pEdgeF * v[ix, jy, kz]

                gRho = flux_muscl(vSurf, rhoB, rhoF)

                rhoFlux[ix, jy, kz, 2] = gRho
            end
        end
    end

    # vertical rho flux
    for kz = 0:nz
        for jy = 1:ny
            for ix = 1:nx
                rhoStratEdgeU =
                    (
                        jac(ix, jy, kz + 1) * rhoStrat[ix, jy, kz] +
                        jac(ix, jy, kz) * rhoStrat[ix, jy, kz+1]
                    ) / (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pEdgeU =
                    (
                        jac(ix, jy, kz + 1) * pStrat[ix, jy, kz] +
                        jac(ix, jy, kz) * pStrat[ix, jy, kz+1]
                    ) / (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                rhoU = rhoTilde[ix, jy, kz+1, 3, 0] + rhoStratEdgeU / pEdgeU
                rhoD = rhoTilde[ix, jy, kz, 3, 1] + rhoStratEdgeU / pEdgeU

                pEdgeU =
                    jac(ix, jy, kz) *
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

function rhopFlux!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz, lz) = grid
    (; var, var0, rhopTilde, pStrat, flux, jac) = cache

    rhopFlux = flux.rhop # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    for kz = 1:nz
        for jy = 1:ny
            for ix = 0:nx
                rhopR = rhopTilde[ix+1, jy, kz, 1, 0]
                rhopL = rhopTilde[ix, jy, kz, 1, 1]

                pEdgeR =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz]
                    )
                uSurf = pEdgeR * u[ix, jy, kz]

                fRhop = flux_muscl(uSurf, rhopL, rhopR)

                rhopFlux[ix, jy, kz, 1] = fRhop

            end
        end
    end

    # meridional rho flux
    for kz = 1:nz
        for jy = 0:ny
            for ix = 1:nx
                rhopF = rhopTilde[ix, jy+1, kz, 2, 0]
                rhopB = rhopTilde[ix, jy, kz, 2, 1]

                pEdgeF =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz]
                    )
                vSurf = pEdgeF * v[ix, jy, kz]

                gRhop = flux_muscl(vSurf, rhopB, rhopF)

                rhopFlux[ix, jy, kz, 2] = gRhop

            end
        end
    end

    # vertical rho flux
    for kz = 0:nz
        for jy = 1:ny
            for ix = 1:nx
                rhopU = rhopTilde[ix, jy, kz+1, 3, 0]
                rhopD = rhopTilde[ix, jy, kz, 3, 1]

                pEdgeU =
                    jac(ix, jy, kz) *
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

function rhouuFlux!(semi) # (rho*u)*u

    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, var0, uTilde, pStrat, flux, jac) = cache

    uFlux = flux.u # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux fRhoU
    for kz = 1:nz
        for jy = 1:ny
            for ix = -1:nx
                uR = uTilde[ix+1, jy, kz, 1, 0]
                uL = uTilde[ix, jy, kz, 1, 1]

                pEdge =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz]
                    )
                pREdge =
                    0.5 * (
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz] +
                        jac(ix + 2, jy, kz) * pStrat[ix+2, jy, kz]
                    )
                uSurf = 0.5 * (pEdge * u[ix, jy, kz] + pREdge * u[ix+1, jy, kz])

                fRhoU = flux_muscl(uSurf, uL, uR)

                uFlux[ix, jy, kz, 1] = fRhoU
            end
        end
    end
end

function rhouvFlux!(semi) # (rho*u)*v

    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, var0, uTilde, pStrat, flux, jac) = cache

    uFlux = flux.u # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux gRhoU
    for kz = 1:nz
        for jy = 0:ny
            for ix = 0:nx
                uF = uTilde[ix, jy+1, kz, 2, 0]
                uB = uTilde[ix, jy, kz, 2, 1]

                pEdge =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz]
                    )
                pREdge =
                    0.5 * (
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz] +
                        jac(ix + 1, jy + 1, kz) * pStrat[ix+1, jy+1, kz]
                    )
                vSurf = 0.5 * (pEdge * v[ix, jy, kz] + pREdge * v[ix+1, jy, kz])

                gRhoU = flux_muscl(vSurf, uB, uF)

                uFlux[ix, jy, kz, 2] = gRhoU
            end
        end
    end
end

function rhouwFlux!(semi) # (rho*u)*w

    (; grid, cache) = semi
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = grid
    (; var, var0, uTilde, pStrat, flux, jac) = cache

    uFlux = flux.u # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux hRhoU
    for kz = 0:nz
        for jy = 1:ny
            for ix = 0:nx
                uU = uTilde[ix, jy, kz+1, 3, 0]
                uD = uTilde[ix, jy, kz, 3, 1]

                pEdge =
                    jac(ix, jy, kz) *
                    jac(ix, jy, kz + 1) *
                    (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                    (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pREdge =
                    jac(ix + 1, jy, kz) *
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


function rhovuFlux!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; var, var0, vTilde, pStrat, flux, jac) = cache

    vFlux = flux.v # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux fRhoV
    for kz = 1:nz
        for jy = 0:ny
            for ix = 0:nx
                vR = vTilde[ix+1, jy, kz, 1, 0]
                vL = vTilde[ix, jy, kz, 1, 1]

                pEdge =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz]
                    )
                pREdge =
                    0.5 * (
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz] +
                        jac(ix + 1, jy + 1, kz) * pStrat[ix+1, jy+1, kz]
                    )
                uSurf = 0.5 * (pEdge * u[ix, jy, kz] + pREdge * u[ix, jy+1, kz])

                fRhoV = flux_muscl(uSurf, vL, vR)

                vFlux[ix, jy, kz, 1] = fRhoV
            end
        end
    end
end

function rhovvFlux!(semi) # (rho*v)*v

    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; var, var0, vTilde, pStrat, flux, jac) = cache

    vFlux = flux.v # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux gRhoV
    for kz = 1:nz
        for jy = -1:ny
            for ix = 1:nx
                vF = vTilde[ix, jy+1, kz, 2, 0]
                vB = vTilde[ix, jy, kz, 2, 1]

                pEdge =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz]
                    )
                pREdge =
                    0.5 * (
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz] +
                        jac(ix, jy + 2, kz) * pStrat[ix, jy+2, kz]
                    )
                vSurf = 0.5 * (pEdge * v[ix, jy, kz] + pREdge * v[ix, jy+1, kz])

                gRhoV = flux_muscl(vSurf, vB, vF)

                vFlux[ix, jy, kz, 2] = gRhoV
            end
        end
    end
end

function rhovwFlux!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; var, var0, vTilde, pStrat, flux, jac) = cache

    vFlux = flux.v # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux hRhoV
    for kz = 0:nz
        for jy = 0:ny
            for ix = 1:nx
                vU = vTilde[ix, jy, kz+1, 3, 0]
                vD = vTilde[ix, jy, kz, 3, 1]

                pEdge =
                    jac(ix, jy, kz) *
                    jac(ix, jy, kz + 1) *
                    (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                    (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pREdge =
                    jac(ix, jy + 1, kz) *
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

function rhowuFlux!(semi) # (rho*w)*u

    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; var, var0, wTilde, pStrat, flux, jac) = cache

    wFlux = flux.w # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux fRhoW
    for kz = 0:nz
        for jy = 1:ny
            for ix = 0:nx
                wR = wTilde[ix+1, jy, kz, 1, 0]
                wL = wTilde[ix, jy, kz, 1, 1]

                pEdge =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix + 1, jy, kz) * pStrat[ix+1, jy, kz]
                    )
                pREdge =
                    0.5 * (
                        jac(ix, jy, kz + 1) * pStrat[ix, jy, kz+1] +
                        jac(ix + 1, jy, kz + 1) * pStrat[ix+1, jy, kz+1]
                    )
                uSurf =
                    (
                        (jac[ix, jy, kz+1] + jac[ix+1, jy, kz+1]) * pEdge * u[ix, jy, kz] +
                        (jac[ix, jy, kz] + jac[ix+1, jy, kz]) * pREdge * u[ix, jy, kz+1]
                    ) / (
                        jac[ix, jy, kz] +
                        jac[ix+1, jy, kz] +
                        jac[ix, jy, kz+1] +
                        jac[ix+1, jy, kz+1]
                    )

                fRhoW = flux_muscl(uSurf, wL, wR)

                wFlux[ix, jy, kz, 1] = fRhoW

            end
        end
    end
end

function rhowvFlux!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; var, var0, wTilde, pStrat, flux, jac) = cache

    wFlux = flux.w # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux gRhoW
    for kz = 0:nz
        for jy = 0:ny
            for ix = 1:nx
                wF = wTilde[ix, jy+1, kz, 2, 0]
                wB = wTilde[ix, jy, kz, 2, 1]

                pEdge =
                    0.5 * (
                        jac(ix, jy, kz) * pStrat[ix, jy, kz] +
                        jac(ix, jy + 1, kz) * pStrat[ix, jy+1, kz]
                    )
                pREdge =
                    0.5 * (
                        jac(ix, jy, kz + 1) * pStrat[ix, jy, kz+1] +
                        jac(ix, jy + 1, kz + 1) * pStrat[ix, jy+1, kz+1]
                    )
                vSurf =
                    (
                        (jac(ix, jy, kz + 1) + jac(ix, jy + 1, kz + 1)) *
                        pEdge *
                        v[ix, jy, kz] +
                        (jac(ix, jy, kz) + jac(ix, jy + 1, kz)) * pREdge * v[ix, jy, kz+1]
                    ) / (
                        jac(ix, jy, kz) +
                        jac(ix, jy + 1, kz) +
                        jac(ix, jy, kz + 1) +
                        jac(ix, jy + 1, kz + 1)
                    )

                gRhoW = flux_muscl(vSurf, wB, wF)

                wFlux[ix, jy, kz, 2] = gRhoW

            end
        end
    end
end

function rhowwFlux!(semi)

    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; var, var0, wTilde, pStrat, flux, jac) = cache

    wFlux = flux.w # mass flux
    u = var0.u # zonal velocity
    v = var0.v # meridional velocity
    w = var0.w # vertical velocity

    # Flux hRhoW
    for kz = -1:nz
        for jy = 1:ny
            for ix = 1:nx
                wU = wTilde[ix, jy, kz+1, 3, 0]
                wD = wTilde[ix, jy, kz, 3, 1]

                pEdge =
                    jac(ix, jy, kz) *
                    jac(ix, jy, kz + 1) *
                    (pStrat[ix, jy, kz] + pStrat[ix, jy, kz+1]) /
                    (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pREdge =
                    jac(ix, jy, kz + 1) *
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

function setHalosOfField!(field, semi)

    (; grid) = semi
    (; nbx, nx, ny, nby) = grid

    for i = 1:nbx
        field[nx+i, :, :] .= field[i, :, :]
        field[-i+1, :, :] .= field[-i+nx+1, :, :]
    end
    for j = 1:nby
        field[:, ny+j, :] .= field[:, j, :]
        field[:, -j+1, :] .= field[:, ny-j+1, :]
    end


end
