#--------- RECONSTRUCTION ---------#
function reconstruct_rho!(var, rhoTilde, grid)

	(; nx, nz, nbx, nbz) = grid

    # rhoTilde is global in PinCFlow
    #rhop = @view var[1, :, :] # density fluctuation
    #rho = rhop + rhoStrat # rhoStrat = background density

	rho = @view var[1, :, :]
    sizeX, sizeZ = size(rho)

    rhoBar_ = zeros(sizeX, sizeZ)
	rhoBar = OffsetArray(rhoBar_, OffsetArrays.Origin(-nbx, -nbz))

    for ix in -nbx:nx+nbx
        for kz in 0:nz+1
            # pStrat is a global variable in PinCFlow
			# TODO: divide by pStrat
            rhoBar[ix, kz] = rho[ix, kz] / 1. # pStrat[ix, kz]
        end
    end
    muscl_reconstruct2D!(rhoBar, nx, nz, rhoTilde)
end

function reconstruct_u!(var, uTilde, grid)

	(; nx, nz, nbx, nbz) = grid

    # uTilde is global in PinCFlow
    rho = @view var[1, :, :] # density fluctuation
    uu = @view var[2, :, :] # zonal velocity
    sizeX, sizeZ = size(uu)

    uBar_ = zeros(sizeX, sizeZ)
	uBar = OffsetArray(uBar_, OffsetArrays.Origin(-nbx, -nbz))

    # rho*u/P for reconstruction
    for ix in -nbx:nx+nbx-1
        for kz in 0:nz+1
            # pStrat is a global variable in PinCFlow
            rhoEdgeR = 0.5 * (rho[ix, kz] + rho[ix + 1, kz])
            # pEdgeR = 0.5 * (pStrat[ix, kz] + pStrat[ix + 1, kz])
            uBar[ix, kz] = uu[ix, kz] * rhoEdgeR / 1. # pEdgeR
        end 
    end
    muscl_reconstruct2D!(uBar, nx, nz, uTilde)
end

function reconstruct_w!(var, wTilde, grid)

	(; nx, nz, nbx, nbz) = grid

    # wTilde is global in PinCFlow
    rhop = @view var[1, :, :] # density fluctuation
    rho = rhop # + rhoStrat # rhoStrat = background density
    ww = @view var[3, :, :] # vertical velocity
    sizeX, sizeZ = size(ww)
    wBar_ = zeros(sizeX, sizeZ)
	wBar = OffsetArray(wBar_, OffsetArrays.Origin(-nbx, -nbz))
    #wBar[:, :] = ww[:, :]

    # for ix in 1:nx 
    #     for kz in 0:nz+1
    #         wBar[ix, kz] = vertWindTFC(var, ix, kz)
    #     end 
    # end
    # TODO: PinCFlow sets halos of wBar
    for ix in -nbx:nx+nbx
        for kz in 0:nz+1
            rhoEdgeU = (jac(ix, kz+1) * rho[ix, kz] + 
            jac(ix, kz) * rho[ix, kz+1]) / (jac(ix, kz) + jac(ix, kz+1))
            #pEdgeU = (jac[ix, kz+1] * pStrat[ix, kz] + 
            #jac[ix, kz] * pStrat[ix, kz+1]) / (jac[ix, kz] + jac[ix, kz+1])
            wBar[ix, kz] = ww[ix, kz] * rhoEdgeU / 1. # pEdgeU
        end
    end
    # muscl_reconstruct2D!(wBar, nx, nz, wTilde)
end

function muscl_reconstruct2D!(var, sizeX, sizeZ, varTilde)

    for k in 0:sizeZ
		orientation = 1
        phiX = var[:, k]
        varTildeX = @view varTilde[:, k, orientation, :] # i, k, orientation, L/R Tilde
        muscle_reconstruct1D_mcvariant!(phiX, sizeX, varTildeX)
    end

    for i in 0:sizeX
        phiZ = var[i, :]
		orientation = 3
        varTildeZ = @view varTilde[i, :, orientation, :] # i, k, orientation, L/R Tilde
        muscle_reconstruct1D_mcvariant!(phiZ, sizeZ, varTildeZ)
    end

end

function muscle_reconstruct1D_mcvariant!(phi, phiSize, phiTilde)

    for i = 0:phiSize
        deltaL = phi[i] - phi[i - 1]
        deltaR = phi[i + 1] - phi[i]

        if (deltaR == 0.)
            phiTilde[i, 2] = phi[i] # 2: phiTilde_R
            phiTilde[i, 1] = phi[i] # 1: phiTilde_L
        else
            if (deltaL == 0.)
                theta = deltaL / deltaR
                s = (2.0 + theta) / 3.0
                sigmaL = max(0., min(2 * theta, s, 2.))

                phiTilde[i, 2] = phi[i] + 0.5 * sigmaL * deltaR
                phiTilde[i, 1] = phi[i]
            else
                theta = deltaL / deltaR
                s = (2.0 + theta) / 3.0
                sigmaL = max(0., min(2 * theta, s, 2.))
                s = (2.0 + 1.0 / theta) / 3.0
                sigmaR = max(0., min(2 / theta, s, 2.))

                phiTilde[i, 2] = phi[i] + 0.5 * sigmaL * deltaR
                phiTilde[i, 1] = phi[i] - 0.5 * sigmaR * deltaL
            end
        end
    end
end

# ------------------------ FLUX CALCULATION ------------------------
function compute_fluxes_PinCFlow!(cache, grid)
	
	rhoFlux!(cache, grid)
	momemtumFlux!(cache, grid)

end

function rhoFlux!(cache, grid)

    (; nx, nz) = grid
    (; u, rhoTilde, flux) = cache

    rhoFlux = @view flux[1, :, :, :] # mass flux
    uu = @view u[2, :, :] # zonal velocity
    ww = @view u[3, :, :] # vertical velocity

    # zonal rho flux 
    orientation = 1
    for kz in 1:nz
        for ix = 0:nx
            # rhoStratEdgeR = 0.5 * (rhoStrat[ix, kz] + rhoStrat[ix + 1, kz])
            # pEdgeR = 0.5 * (pStrat[ix, kz] + pStrat[ix + 1, kz])
            rhoR = rhoTilde[ix + 1, kz, orientation, 1] 
                # + rhoStratEdgeR / pEdgeR
            rhoL = rhoTilde[ix, kz, orientation, 2] 
                # + rhoStratEdgeR / pEdgeR

            # pEdgeR = 0.5 * (jac(ix, kz) * pStrat[ix, kz] 
            #     + jac(ix + 1, kz) * pStrat[ix + 1, kz])
            pEdgeR = 1.
            uSurf = pEdgeR * uu[ix, kz]

            fRho = flux_muscl(uSurf, rhoL, rhoR)

            rhoFlux[ix, kz, orientation] = fRho
        end
    end

    # vertical rho flux
    orientation = 3
    for kz in 0:nz
        for ix in 1:nx
            # rhoStratEdgeU = (jac(ix, kz+1) * rhoStrat[ix, kz] + 
            #     jac(ix, kz) * rhoStrat(ix, kz+1)) / 
            #     (jac(ix, kz) + jac(ix, kz+1))
            # pEdgeU = (jac(ix, kz+1) * pStrat[ix, kz] + 
            #     jac(ix, kz) * pStrat[ix, kz + 1]) / 
            #     (jac(ix, kz) + jac(ix, kz+1))
            rhoU = rhoTilde[ix, kz + 1, orientation, 1] 
                # + rhoStratEdgeU / pEdgeU
            rhoD = rhoTilde[ix, kz, orientation, 2] 
                # + rhoStratEdgeU / pEdgeU
            # pEdgeU = jac(ix, kz) * jac(ix, kz + 1) * (pStrat[ix, kz] + 
            #     pStrat[ix, kz + 1]) / (jac(ix, kz) + jac(ix, kz + 1))
            pEdgeU = 1.
            wSurf = pEdgeU * ww[ix, kz]

            hRho = flux_muscl(wSurf, rhoD, rhoU)

            rhoFlux[ix, kz, orientation] = hRho
        end
    end
         
end

function momemtumFlux!(cache, grid)

    zonalMomentumFlux!(cache, grid) # (rho*u)*u, (rho*u)*w
    # TODO meridionalMomentumFlux(var, flux)
    verticalMomentumFlux!(cache, grid) # (rho*w)*u, (rho*w)*w

end

function zonalMomentumFlux!(cache, grid)

    rhouuFlux!(cache, grid) # calculate (rho*u)*u
    # TODO rhouvFlux(var, flux) # calculate (rho*u)*v
    rhouwFlux!(cache, grid) # calculate (rho*u)*w

end

function verticalMomentumFlux!(cache, grid)

    rhowuFlux!(cache, grid) # calculate (rho*w)*u
    # TODO rhowvFlux(var, flux) # calculate (rho*w)*v
    rhowwFlux!(cache, grid) # calculate (rho*w)*w

end

function rhouuFlux!(cache, grid) # (rho*u)*u

    (; nx, nz) = grid
    (; u, uTilde, flux) = cache

    orientation = 1
    uFlux = @view flux[2, :, :, orientation] # zonal momentum flux
    uu = @view u[2, :, :] # zonal velocity

    # Flux fRhoU
    for kz in 0:nz
        for ix in 0:nx
            uR = uTilde[ix + 1, kz, orientation, 1]
            uL = uTilde[ix, kz, orientation, 2]

            # pEdgeR = 0.5 * (jac(ix, kz) * pStrat[ix, kz] + 
            #     jac(ix + 1, kz) * pStrat[ix + 1, kz])
            pEdgeR = 1.
            # pREdgeR = 0.5 * (jac(ix + 1, kz) * pStrat[ix + 1, kz] + 
            #     jac(ix + 2, kz) * pStrat[ix + 2, kz])
            pREdgeR = 1.
            uSurf = 0.5 * (pEdgeR * uu[ix, kz] + pREdgeR * uu[ix + 1, kz])

            fRhoU = flux_muscl(uSurf, uL, uR)

            uFlux[ix, kz] = fRhoU
        end
    end
end

function rhouwFlux!(cache, grid) # (rho*u)*w

    (; nx, nz) = grid
    (; u, uTilde, flux) = cache

    orientation = 3
    uFlux = @view flux[2, :, :, orientation] # zonal momentum flux
    ww = @view u[3, :, :] # vertical velocity

    # Flux hRhoU
    for kz in 0:nz
        for ix in 0:nx

            uU = uTilde[ix, kz + 1, orientation, 1]
            uD = uTilde[ix, kz, orientation, 2]

            # pEdge = jac(ix, kz) * jac(ix, kz + 1) * (pStrat[ix, kz] + 
            #     pStrat[ix, kz + 1]) / (jac(ix, kz) + jac(ix, kz + 1))
            pEdge = 1.
            # pREdge = jac(ix + 1, kz) * jac(ix + 1, kz + 1) / (jac(ix + 1, kz) 
            #     + jac(ix + 1, kz + 1))
            pREdge = 1.
            wSurf = 0.5 * (pEdge * ww[ix, kz] + pREdge * ww[ix + 1, kz])

            hRhoU = flux_muscl(wSurf, uD, uU)

            uFlux[ix, kz] = hRhoU
        end 
    end
end

function rhowuFlux!(cache, grid) # (rho*w)*u

    (; nx, nz) = grid
    (; u, wTilde, flux) = cache

    orientation = 1
    wFlux = @view flux[3, :, :, orientation] # vertical momentum flux
    uu = @view u[2, :, :] # zonal velocity

    # Flux fRhoW
    for kz in 0:nz
        for ix in 0:nx
            wR = wTilde[ix + 1, kz, orientation, 1]
            wL = wTilde[ix, kz, orientation, 2]

            # pEdge = 0.5 * (jac(ix, kz) * pStrat[ix, kz] + 
            #     jac(ix + 1, kz) * pStrat[ix + 1, kz])
            pEdge = 1.
            # pUEdge = 0.5 * (jac(ix, kz + 1) * pStrat[ix, kz + 1] + 
            #     jac(ix + 1, kz + 1) * pStrat[ix + 1, kz + 1])
            pUEdge = 1.
            uSurf = ((jac(ix, kz + 1) + jac(ix + 1, kz + 1)) * pEdge 
                * uu[ix, kz] + (jac(ix, kz) + jac(ix + 1, kz)) * pUEdge 
                * uu[ix, kz + 1]) / (jac(ix, kz) + jac(ix + 1, kz) 
                + jac(ix, kz + 1) + jac(ix + 1, kz + 1))
                
            fRhoW = flux_muscl(uSurf, wL, wR)

            wFlux[ix, kz] = fRhoW
        end
    end
end

function rhowwFlux!(cache, grid) # (rho*w)*w

    (; nx, nz) = grid
    (; u, wTilde, flux) = cache

    orientation = 3
    wFlux = @view flux[3, :, :, orientation] # vertical momentum flux
    ww = @view u[3, :, :] # vertical velocity

    # Flux hRhoW
    for kz in 0:nz
        for ix in 0:nx

            # println("(ix, kz) = ", (ix, kz))
            wU = wTilde[ix, kz + 1, orientation, 1]
            wD = wTilde[ix, kz, orientation, 2]

            # TODO: inconsistent naming of pEdge and pUEdge
            # pEdge = jac(ix, kz) * jac(ix, kz + 1) * (pStrat[ix, kz] + 
            #     pStrat[ix, kz + 1]) / (jac(ix, kz) + jac(ix, kz + 1))
            pEdge = 1.
            # pUEdge = jac(ix, kz + 1) * jac(ix, kz + 2) * (pStrat[ix, kz + 1] + 
            #     pStrat[ix, kz + 2]) / (jac(ix, kz + 1) + jac(ix, kz + 2))
            pUEdge = 1.
            wSurf = 0.5 * (pEdge * ww[ix, kz] + pUEdge * ww[ix, kz + 1])

            hRhoW = flux_muscl(wSurf, wD, wU)

            wFlux[ix, kz] = hRhoW
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

function jac(ix, kz)
    # Jacobian
    # TODO: include topography
    return jac = 1.

    # jac = (lz[2] - topography_surface[ix]) / (lz[2] * 
    #   zTildeS[kz] - zTildeS[kz + 1]) / dz
end
