using SimpleUnPack
using OffsetArrays

function make_grid(; nx, ny, nz, nbx, nby, nbz, xmin, xmax, ymin, ymax, zmin, zmax, lRef)
    topography_surface = OffsetArray(zeros(nx + 1 + 2*nbx, ny + 1 + 2*nby), -nbx:nx+nbx, -nby:ny+nby)
    zTildeTFC = OffsetArray(zeros(nx + 1 + 2*nbx, ny + 1 + 2*nby, nz + 1 + 2*nbz), -nbx:nx+nbx, -nby:ny+nby, -nbz:nz+nbz)
    zTFC = OffsetArray(zeros(nx + 1 + 2*nbx, ny + 1 + 2*nby, nz + 1 + 2*nbz), -nbx:nx+nbx, -nby:ny+nby, -nbz:nz+nbz)
    zTildeS = OffsetArray(zeros(nz + 1 + 2*nbz), -nbz:nz+nbz)
    zS = OffsetArray(zeros(nz + 1 + 2*nbz), -nbz:nz+nbz)

    lx_dim = OffsetArray([xmin, xmax], 0:1)
    ly_dim = OffsetArray([ymin, ymax], 0:1)
    lz_dim = OffsetArray([zmin, zmax], 0:1)

    lx = lx_dim ./ lRef
    ly = ly_dim ./ lRef
    lz = lz_dim ./ lRef

    dx = (lx[1] - lx[0])
    dy = (ly[1] - ly[0])
    dz = (lz[1] - lz[0])


    sizeX = nx
    sizeY = ny
    sizeZ = nz
    x = OffsetArray(zeros(sizeX + 1 + 2nbx), -nbx:sizeX+nbx)
    y = OffsetArray(zeros(sizeY + 1 + 2nby), -nby:sizeY+nby)
    z = OffsetArray(zeros(sizeZ + 1 + 2nbz), -nbz:sizeZ+nbz)

    for i in -nbx:sizeX+nbx
        x[i] = lx[0] + (i - 1) * dx + dx / 2.0
    end

    for j in -nby:sizeY+nby
        y[j] = ly[0] + (j - 1) * dy + dy / 2.0
    end

    for k in -nbz:sizeZ+nbz
        z[k] = lz[0] + (k - 1) * dz + dz / 2.0
    end

    grid = (; nx, ny, nz,
             nbx, nby, nbz, topography_surface, zTildeTFC, zTFC, zTildeS, zS,
             lx, ly, lz, dx, dy, dz, x, y, z)
    return grid
end

nx = ny = nz = 100
nbx = nby = nbz = 3

# Constants
gamma = 1.4
gamma_1 = gamma - 1.0
kappa = (gamma - 1.0) / gamma
kappaInv = 1.0 / kappa
gammaInv = 1.0 / gamma

# Reference quantities
Rsp = 287.0
g = 9.81
rhoRef = 1.184
pRef = 101325.0
aRef = sqrt(pRef / rhoRef)
uRef = aRef
global lRef = pRef / rhoRef / g
tRef = lRef / aRef
thetaRef = aRef ^ 2. / Rsp
Ma = uRef / aRef
Fr = uRef / sqrt(g * lRef)
kappa = (gamma - 1.) / gamma
sig = Ma ^ 2 / Fr ^ 2

FRef = rhoRef * uRef^2 / lRef

press0_dim = 1.0e5

# isothermal atmosphere

Temp0_dim = 300.0
T0 = Temp0_dim / thetaRef
N2 = Ma ^ 2 / Fr ^ 4 * kappa / T0
NN = sqrt(N2)

mu_viscous_dim = 0.0
ReInv = mu_viscous_dim / (uRef * lRef)

# Reynolds number
ReInv = 0.0
if false # TODO - Figure it out!
ReInv = mu_viscous_dim / (uRef * lRef)
end

Re = if ReInv < 1.0e-20
1.0e20
else
1.0 / ReInv
end

xmin = ymin = zmin = 0.0
xmax = ymax = zmax = 1.0
grid = make_grid(nx=nx, ny=ny, nz=nz, nbx=nbx, nby=nby, nbz=nbz, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, zmin = zmin, zmax = zmax, lRef = lRef)

equations = (; gamma, gamma_1, kappaInv, gammaInv, Rsp, g, rhoRef, pRef, aRef, uRef, lRef, tRef, thetaRef, Ma, Fr, kappa, sig, press0_dim, Temp0_dim, T0, N2, NN, mu_viscous_dim, ReInv, Re)

semi = (; grid, equations)

function initialize_atmosphere!(semi)

    (; grid, equations) = semi
    (; nx, ny, nz,
             nbx, nby, nbz, topography_surface, zTildeTFC, zTFC, zTildeS, zS,
             lx, ly, lz, dx, dy, dz, x, y, z) = grid

    (; gamma, gamma_1, kappa, kappaInv, gammaInv, Rsp, g, rhoRef, pRef, aRef, uRef, lRef, tRef, thetaRef, Ma, Fr, kappa, sig, press0_dim, Temp0_dim, T0, N2, NN, mu_viscous_dim, ReInv, Re) = equations

    setup_topography!(topography_surface, zTFC, nx, ny, nz, lz)

    # scaled background flow
    backgroundFlow_dim = 10.0
    backgroundFlow = backgroundFlow_dim / uRef

    # nondimensional gravitational constant
    g_ndim = g / (uRef ^ 2 / lRef)
    p0 = press0_dim / pRef
    pStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                        -nbx:nx+nbx, -nby:ny+nby, -1:nz+2)

    thetaStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                            -nbx:nx+nbx, -nby:ny+nby, -1:nz+2)

    rhoStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                        -nbx:nx+nbx, -nby:ny+nby, -1:nz+2)

    bvsStrat = OffsetArray(zeros(Float64, nx + 2nbx + 1, ny + 2nby + 1, nz + 4),
                        -nbx:nx+nbx, -nby:ny+nby, -1:nz+2)



    for ix in -nbx:nx + nbx, jy in -nby:ny + nby
        for kz in -1:nz + 2
            zTFC[ix, jy, kz]
            pStrat[ix, jy, kz] = p0 * exp(- sig * zTFC[ix, jy, kz] / gamma / T0)
            thetaStrat[ix, jy, kz] = T0 * exp(kappa * sig / T0 * zTFC[ix, jy, kz])
            rhoStrat[ix, jy, kz] = pStrat[ix, jy, kz] / thetaStrat[ix, jy, kz]
        end
    end

    for ix in -nbx:nx + nbx, jy in -nby:ny + nby
        bvsStrat[ix, jy, -1] = g_ndim / thetaStrat[ix, jy, 0] / jac(ix, jy, 0, lz, grid) *
            (thetaStrat[ix, jy, 1] - thetaStrat[ix, jy, 0]) / dz
        for kz in 1:nz
            bvsStrat[ix, jy, kz] = g_ndim / thetaStrat[ix, jy, kz] / jac(ix, jy, kz, lz, grid) *
                0.5 * (thetaStrat[ix, jy, kz + 1] - thetaStrat[ix, jy, kz - 1]) / dz
        end
        bvsStrat[ix, jy, nz + 1] = g_ndim / thetaStrat[ix, jy, nz + 1] / jac(ix, jy, nz + 1, lz, grid) *
            (thetaStrat[ix, jy, nz + 1] - thetaStrat[ix, jy, nz]) / dz
        bvsStrat[ix, jy, nz + 2] = bvsStrat[ix, jy, nz + 1]
    end

end

function setup_topography!(topography_surface, zTFC, nx, ny, nz, lz)

    # (; xc, zc, xf, zf, nx, nz) = grid

    (; grid, equations) = semi
    (; nx, ny, nz,  nbx, nby, nbz, topography_surface, zTildeTFC, zTFC, zTildeS, zS,
             lx, ly, lz, dx, dy, dz, x, y, z) = grid

    (; gamma, gamma_1, kappa, kappaInv, gammaInv, Rsp, g, rhoRef, pRef, aRef, uRef, lRef, tRef, thetaRef, Ma, Fr, kappa, sig, press0_dim, Temp0_dim, T0, N2, NN, mu_viscous_dim, ReInv, Re) = equations

    if lz[0] != 0.0
        @assert false "Error in setup_topography: lz(0) must be zero for & &TFC!"
    end

    mountainHeight_dim = mountainWidth_dim = 1.0 # TODO - PLEASE!!
    mountainHeight = mountainHeight_dim / lRef
    mountainWidth = mountainWidth_dim / lRef
    mountainWavenumber = pi / mountainWidth

    x_center = 0.5 * (lx[1] + lx[0])
    y_center = 0.5 * (ly[1] + ly[0])
    mountain_case = 2

    if mountain_case != 0
      topography_surface .= 0.0
      for jy = 1:ny
        for ix = 1:nx
          topography_surface[ix, jy] = mountainHeight / (1.0 + (x[ix] - x_center) ^ 2.0 / mountainWidth ^ 2.0)
        end
      end
    else
      topography_surface = topography_surface / lRef
    end

    # TODO - check Halos
    # setHalosOfField2D(topography_surface)

    # Compute the stretched vertical grid.
    for kz = - nbz:nz + nbz
      zTildeS[kz] = map(z[kz] + 0.5 * dz, lz)
    end
    for kz = - nbz + 1: nz + nbz
      zS[kz] = 0.5 * (zTildeS[kz] + zTildeS[kz - 1])
    end
    zS[- nbz] = zTildeS[- nbz] - 0.5 * (zTildeS[nbz + 1] - zTildeS[nbz])

    # Compute the physical layers.
    for kz = - nbz:nz + nbz
      zTFC[:, :, kz] .= (lz[1] .- topography_surface) / lz[1] * zS[kz] .+ topography_surface
    end
end

function map(level, lz)
    stretch_exponent = 2
    # Vertical grid stretching.
    if level < 0
        return -lz[1] * (-level / lz[1])^stretch_exponent
    elseif level > lz[1]
        return 2 * lz[1] - lz[1] * ((2 * lz[1] - level) / lz[1])^stretch_exponent
    else
        return lz[1] * (level / lz[1])^stretch_exponent
    end
end

function jac(i, j, k, lz, grid)
    # Jacobian.
    (; topography_surface, zTildeS, dz) = grid
    return (lz[1] - topography_surface[i, j]) / lz[1] * (zTildeS[k] - zTildeS[k - 1]) / dz
end

function met(i, j, k, mu, nu)
    # Metric tensor.

    if (mu == 1 && nu == 3) || (mu == 3 && nu == 1)
        return (topography_surface[i + 1, j] - topography_surface[i - 1, j]) /
               (2.0 * dx) * (zS[k] - lz[1]) / (lz[1] - topography_surface[i, j]) *
               dz / (zTildeS[k] - zTildeS[k - 1])
    elseif (mu == 2 && nu == 3) || (mu == 3 && nu == 2)
        return (topography_surface[i, j + 1] - topography_surface[i, j - 1]) /
               (2.0 * dy) * (zS[k] - lz[1]) / (lz[1] - topography_surface[i, j]) *
               dz / (zTildeS[k] - zTildeS[k - 1])
    elseif mu == 3 && nu == 3
        return ((lz[1] / (lz[1] - topography_surface[i, j]))^2.0 +
                ((zS[k] - lz[1]) / (lz[1] - topography_surface[i, j]))^2.0 *
                (((topography_surface[i + 1, j] - topography_surface[i - 1, j]) / (2.0 * dx))^2.0 +
                 ((topography_surface[i, j + 1] - topography_surface[i, j - 1]) / (2.0 * dy))^2.0)) *
               (dz / (zTildeS[k] - zTildeS[k - 1]))^2.0
    else
        @assert false "UNDEFINED CASE!!!"
    end
end

function vertWind(i, j, k, var)
    # Transformation of the vertical wind.

    uEdgeR = var.u[i, j, k]
    uUEdgeR = var.u[i, j, k + 1]
    uEdgeL = var.u[i - 1, j, k]
    uUEdgeL = var.u[i - 1, j, k + 1]
    vEdgeF = var.v[i, j, k]
    vUEdgeF = var.v[i, j, k + 1]
    vEdgeB = var.v[i, j - 1, k]
    vUEdgeB = var.v[i, j - 1, k + 1]
    wEdgeU = var.w[i, j, k]

    return trafo(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF,
                 vUEdgeF, vEdgeB, vUEdgeB, wEdgeU, "car")
end

function trafo(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB, vUEdgeB, wEdgeU, wind)
    # Assuming jac and met are defined elsewhere
    # Define variables as in the original code
    jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) + jac(i, j, k + 1))
    uC = 0.5 * (uEdgeR + uEdgeL)
    uU = 0.5 * (uUEdgeR + uUEdgeL)
    vC = 0.5 * (vEdgeF + vEdgeB)
    vU = 0.5 * (vUEdgeF + vUEdgeB)

    if wind == "car"
        trafo = jacEdgeU * (-(jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC + met(i, j, k, 2, 3) * vC) +
            jac(i, j, k) * (met(i, j, k + 1, 1, 3) * uU + met(i, j, k + 1, 2, 3) * vU)) / (jac(i, j, k) + jac(i, j, k + 1)) + wEdgeU)
    elseif wind == "tfc"
        trafo = (jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC + met(i, j, k, 2, 3) * vC) +
            jac(i, j, k) * (met(i, j, k + 1, 1, 3) * uU + met(i, j, k + 1, 2, 3) * vU)) / (jac(i, j, k) + jac(i, j, k + 1)) + wEdgeU / jacEdgeU
    end
    return trafo
end

function stressTensTFC(i, j, k, mu, nu, var)
    # Assuming jac, met, and vertWindTFC are functions defined elsewhere

    # Define variables as in the original code
    jacEdgeR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
    jacEdgeL = 0.5 * (jac(i, j, k) + jac(i - 1, j, k))
    jacEdgeF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k))
    jacEdgeB = 0.5 * (jac(i, j, k) + jac(i, j - 1, k))
    jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) + jac(i, j, k + 1))
    jacEdgeD = 2.0 * jac(i, j, k) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k - 1))

    # Accessing array elements in var
    uF = 0.5 * (var.u[i, j + 1, k] + var.u[i - 1, j + 1, k])
    uB = 0.5 * (var.u[i, j - 1, k] + var.u[i - 1, j - 1, k])
    uU = 0.5 * (var.u[i, j, k + 1] + var.u[i - 1, j, k + 1])
    uD = 0.5 * (var.u[i, j, k - 1] + var.u[i - 1, j, k - 1])

    vR = 0.5 * (var.v[i + 1, j, k] + var.v[i + 1, j - 1, k])
    vL = 0.5 * (var.v[i - 1, j, k] + var.v[i - 1, j - 1, k])
    vU = 0.5 * (var.v[i, j, k + 1] + var.v[i, j - 1, k + 1])
    vD = 0.5 * (var.v[i, j, k - 1] + var.v[i, j - 1, k - 1])

    wR = 0.5 * (vertWindTFC(i + 1, j, k, var) + vertWindTFC(i + 1, j, k - 1, var))
    wL = 0.5 * (vertWindTFC(i - 1, j, k, var) + vertWindTFC(i - 1, j, k - 1, var))
    wF = 0.5 * (vertWindTFC(i, j + 1, k, var) + vertWindTFC(i, j + 1, k - 1, var))
    wB = 0.5 * (vertWindTFC(i, j - 1, k, var) + vertWindTFC(i, j - 1, k - 1, var))

    # Conditional logic for stress tensor calculation
    if mu == 1 && nu == 1
        stressTensTFC = 2.0 * (var.u[i, j, k] - var.u[i - 1, j, k]) / dx +
            met(i, j, k, 1, 3) * (uU - uD) / dz - 2.0 / 3.0 * ((jacEdgeR * var.u[i, j, k] - jacEdgeL * var.u[i - 1, j, k]) / dx +
            (jacEdgeF * var.v[i, j, k] - jacEdgeB * var.v[i, j - 1, k]) / dy +
            (jacEdgeU * var.w[i, j, k] - jacEdgeD * var.w[i, j, k - 1]) / dz) / jac(i, j, k)
    elseif (mu == 1 && nu == 2) || (mu == 2 && nu == 1)
        stressTensTFC = 0.5 * (uF - uB) / dy + 0.5 * met(i, j, k, 2, 3) * (uU - uD) / dz + 0.5 * (vR - vL) / dx +
            0.5 * met(i, j, k, 1, 3) * (vU - vD) / dz
    elseif (mu == 1 && nu == 3) || (mu == 3 && nu == 1)
        stressTensTFC = 0.5 * (uU - uD) / dz / jac(i, j, k) + 0.5 * (wR - wL) / dx + met(i, j, k, 1, 3) * (vertWindTFC(i, j, k, var) - vertWindTFC(i, j, k - 1, var)) / dz
    elseif mu == 2 && nu == 2
        stressTensTFC = 2.0 * (var.v[i, j, k] - var.v[i, j - 1, k]) / dy +
            met(i, j, k, 2, 3) * (vU - vD) / dz - 2.0 / 3.0 * ((jacEdgeR * var.u[i, j, k] - jacEdgeL * var.u[i - 1, j, k]) / dx +
            (jacEdgeF * var.v[i, j, k] - jacEdgeB * var.v[i, j - 1, k]) / dy +
            (jacEdgeU * var.w[i, j, k] - jacEdgeD * var.w[i, j, k - 1]) / dz) / jac(i, j, k)
    elseif (mu == 2 && nu == 3) || (mu == 3 && nu == 2)
        stressTensTFC = 0.5 * (vU - vD) / dz / jac(i, j, k) + 0.5 * (wF - wB) / dy + met(i, j, k, 2, 3) * (vertWindTFC(i, j, k, var) - vertWindTFC(i, j, k - 1, var)) / dz
    elseif mu == 3 && nu == 3
        stressTensTFC = 2.0 * (vertWindTFC(i, j, k, var) - vertWindTFC(i, j, k - 1, var)) / dz / jac(i, j, k) -
            2.0 / 3.0 * ((jacEdgeR * var.u[i, j, k] - jacEdgeL * var.u[i - 1, j, k]) / dx +
            (jacEdgeF * var.v[i, j, k] - jacEdgeB * var.v[i, j - 1, k]) / dy +
            (jacEdgeU * var.w[i, j, k] - jacEdgeD * var.w[i, j, k - 1]) / dz) / jac(i, j, k)
    end

    return stressTensTFC
end

initialize_atmosphere!(semi)
