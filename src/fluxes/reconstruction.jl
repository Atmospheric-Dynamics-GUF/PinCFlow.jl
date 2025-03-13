# TODO: check start and end of loops
function reconstruction!(model)
    @trixi_timeit timer() "Reconstruction" begin
    #! format: noindent
    reconstruct_rho!(model)
    reconstruct_rhop!(model)
    reconstruct_u!(model)
    reconstruct_v!(model)
    reconstruct_w!(model)
    end # timer
end

function reconstruct_rho!(model)
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = model.domain

    rho = model.variables.prognostic_fields.rho
    pstrat = model.atmosphere.pstrattfc
    sizeX, sizeY, sizeZ = size(rho)

    rhoBar_ = zeros(sizeX, sizeY, sizeZ)
    rhoBar = OffsetArray(rhoBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    for ix in (-nbx):(nx+nbx)
        for jy in (-nby):(ny+nby)
            for kz in 0:(nz+1)
                # TODO: error message for pStrat
                rhoBar[ix, jy, kz] = rho[ix, jy, kz] / pstrat[ix, jy, kz]
            end
        end
    end
    muscl_reconstruct3D!(rhoBar, nxx, nyy, nzz, model.variables.reconstructed.rho)
end

function reconstruct_rhop!(model)
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = model.domain
    rhop = model.variables.prognostic_fields.rhop
    pstrat = model.atmosphere.pstrattfc
    sizeX, sizeY, sizeZ = size(rhop)

    rhopBar_ = zeros(sizeX, sizeY, sizeZ)
    rhopBar = OffsetArray(rhopBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    for ix in (-nbx):(nx+nbx)
        for jy in (-nby):(ny+nby)
            for kz in 0:(nz+1)
                rhopBar[ix, jy, kz] = rhop[ix, jy, kz] / pstrat[ix, jy, kz]
            end
        end
    end
    muscl_reconstruct3D!(rhopBar, nxx, nyy, nzz, model.variables.reconstructed.rhop)
    # TODO: keep limiterType1 ?
end

function reconstruct_u!(model)
    # Compute rho*u/P for reconstruction
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = model.domain
    (; u, rho) = model.variables.prognostic_fields
    pstrat = model.atmosphere.pstrattfc
    rhostrat = model.atmosphere.rhostrattfc
    sizeX, sizeY, sizeZ = size(u)

    uBar_ = zeros(sizeX, sizeY, sizeZ)
    uBar = OffsetArray(uBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    # rho*u/P for reconstruction
    for ix in (-nbx):(nx+nbx-1)
        for jy in (-nby):(ny+nby)
            for kz in 0:(nz+1)
                rhoEdge = 0.5 * (rho[ix, jy, kz] +
                                 rho[ix+1, jy, kz] +
                                 rhostrat[ix, jy, kz] +
                                 rhostrat[ix+1, jy, kz])
                pEdge = 0.5 * (pstrat[ix, jy, kz] + pstrat[ix+1, jy, kz])
                uBar[ix, jy, kz] = u[ix, jy, kz] * rhoEdge / pEdge
            end
        end
    end

    muscl_reconstruct3D!(uBar, nxx, nyy, nzz, model.variables.reconstructed.u)
end

function reconstruct_v!(model)

    # Compute rho*v/P for reconstruction
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = model.domain
    (; u, v, rho) = model.variables.prognostic_fields
    pstrat = model.atmosphere.pstrattfc
    rhostrat = model.atmosphere.rhostrattfc

    sizeX, sizeY, sizeZ = size(u)

    vBar_ = zeros(sizeX, sizeY, sizeZ)
    vBar = OffsetArray(vBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    # rho*v/P for reconstruction
    for ix in (-nbx):(nx+nbx)
        for jy in (-nby):(ny+nby-1)
            for kz in 0:(nz+1)
                rhoEdge = 0.5 * (rho[ix, jy, kz] +
                                 rhostrat[ix, jy, kz] +
                                 rho[ix, jy+1, kz] +
                                 rhostrat[ix, jy+1, kz])
                pEdge = 0.5 * (pstrat[ix, jy, kz] + pstrat[ix, jy+1, kz])
                vBar[ix, jy, kz] = v[ix, jy, kz] * rhoEdge / pEdge
            end
        end
    end

    muscl_reconstruct3D!(vBar, nxx, nyy, nzz, model.variables.reconstructed.v)
end

function reconstruct_w!(model)

    # Compute rho*v/P for reconstruction
    (; nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz) = model.domain
    (; u, v, w, rho) = model.variables.prognostic_fields
    jac = model.grid.jac
    pstrat = model.atmosphere.pstrattfc
    rhostrat = model.atmosphere.rhostrattfc
    sizeX, sizeY, sizeZ = size(w)

    wBar_ = zeros(sizeX, sizeY, sizeZ)
    wBar = OffsetArray(wBar_, OffsetArrays.Origin(-nbx, -nby, -nbz))

    wBar[:, :, 0:(nz+1)] .= w[:, :, 0:(nz+1)]
    # TODO: vertWind doesn't seem to work!!!
    for ix in 1:nx
        for jy in 1:ny
            for kz in 0:(nz+1)
                wBar[ix, jy, kz] = vertWind(ix, jy, kz, model)
            end
        end
    end

    setHalosOfField!(wBar, model)
    for ix in (-nbx):(nx+nbx)
        for jy in (-nby):(ny+nby)
            for kz in 0:(nz+1)
                rhoEdge = (jac(ix, jy, kz + 1) * (rho[ix, jy, kz] + rhostrat[ix, jy, kz]) +
                           jac(ix, jy, kz) *
                           (rho[ix, jy, kz+1] + rhostrat[ix, jy, kz+1])) /
                          (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                pEdge = (jac(ix, jy, kz + 1) * pstrat[ix, jy, kz] +
                         jac(ix, jy, kz) * pstrat[ix, jy, kz+1]) /
                        (jac(ix, jy, kz) + jac(ix, jy, kz + 1))
                wBar[ix, jy, kz] = wBar[ix, jy, kz] * rhoEdge / pEdge
            end
        end
    end
    muscl_reconstruct3D!(wBar, nxx, nyy, nzz, model.variables.reconstructed.w; debug=1)
end

function muscl_reconstruct3D!(varBar, sizeX, sizeY, sizeZ, varTilde; debug=0)

    # TODO: changed indexing. Check if output is correct!

    phiX = zeros(sizeX)
    phiY = zeros(sizeY)
    phiZ = zeros(sizeZ)

    phiTildeX = OffsetArray(zeros(sizeX, 2), 1:sizeX, 0:1)
    phiTildeY = OffsetArray(zeros(sizeY, 2), 1:sizeY, 0:1)
    phiTildeZ = OffsetArray(zeros(sizeZ, 2), 1:sizeZ, 0:1)

    for kz in 2:(sizeZ-1)
        for jy in 2:(sizeY-1)
            phiX[:] .= varBar.parent[:, jy, kz]
            muscle_reconstruct1D_mcvariant!(phiX, sizeX, phiTildeX)
            varTilde.parent[:, jy, kz, 1, :] .= phiTildeX.parent
        end
    end

    for kz in 2:(sizeZ-1)
        for ix in 2:(sizeX-1)
            phiY .= varBar.parent[ix, :, kz]
            muscle_reconstruct1D_mcvariant!(phiY, sizeY, phiTildeY)
            varTilde.parent[ix, :, kz, 2, :] .= phiTildeY.parent
        end
    end

    for jy in 2:(sizeY-1)
        for ix in 2:(sizeX-1)
            phiZ .= varBar.parent[ix, jy, :]
            muscle_reconstruct1D_mcvariant!(phiZ, sizeZ, phiTildeZ)
            varTilde.parent[ix, jy, :, 3, :] .= phiTildeZ.parent
        end
    end
end

function muscle_reconstruct1D_mcvariant!(phi, phiSize, phiTilde)
    for i in 2:(phiSize-1)
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


function setHalosOfField!(field, model)
    (; nbx, nx, ny, nby) = model.domain

    for i in 1:nbx
        field[nx+i, :, :] .= field[i, :, :]
        field[-i+1, :, :] .= field[-i+nx+1, :, :]
    end
    for j in 1:nby
        field[:, ny+j, :] .= field[:, j, :]
        field[:, -j+1, :] .= field[:, ny-j+1, :]
    end
end
