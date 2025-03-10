struct Time{I<:Integer,VOF<:Vector{<:AbstractFloat}}
    nstages::I
    alphark::VOF
    betark::VOF
    stepfrac::VOF
end

function Time()
    nstages = 3
    betaRK = [1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0]
    alphaRK = [0.0, -5.0 / 9.0, -153.0 / 128.0]
    stepfrac = [1 / 3, 5 / 12, 1 / 4]
    return Time(nstages, alphaRK, betaRK, stepfrac)
end

struct Model
    parameters::Parameters
    constants::Constants
    time::Time
    domain::Domain
    grid::Grid
    atmosphere::Atmosphere
    operator::PoissonOperator
    variables::Variables
    fluxes::Fluxes
end

function Model(p::Parameters)
    c = Constants(p)
    g = Grid(p, c.lref)
    d = Domain(p)
    t = Time()
    a = Atmosphere(p, c)
    o = PoissonOperator(p)
    fl = Fluxes(p.domain)
    v = Variables(p)
    return Model(p, c, t, d, g, a, o, v, fl)
end


function initialize!(model::Model)
    initialize_atmosphere!(model)

end

function setup_topography!(model::Model)
    @trixi_timeit timer() "Setup topography" begin
    #! format: noindent

    # (; xc, zc, xf, zf, nx, nz) = grid

    (; grid, variables, domain) = model
    (; nx, ny, nz) = domain
    (; nbz) = model.parameters.domain
    (;
        topography_surface, ztfc, ztildes, zs, lx, ly, lz, dx, dy, dz, x, y, z,) = grid
    c = model.constants
    if lz[0] != 0.0
        @assert false "Error in setup_topography: lz(0) must be zero for & &TFC!"
    end

    mountainHeight = model.parameters.topography.mountainheight_dim / c.lref
    mountainWidth = model.parameters.topography.mountainwidth_dim / c.lref
    mountainWavenumber = pi / mountainWidth

    x_center = 0.5 * (lx[1] + lx[0])
    y_center = 0.5 * (ly[1] + ly[0])
    mountain_case = model.parameters.topography.mountain_case
    if mountain_case != 0
        topography_surface .= 0.0
        for jy in 1:ny
            for ix in 1:nx
                topography_surface[ix, jy] = mountainHeight / (1.0 +
                                                               (x[ix] - x_center)^2.0 /
                                                               mountainWidth^2.0)
            end
        end
    else
        topography_surface = topography_surface / c.lref
    end

    # TODO - check Halos
    # setHalosOfField2D(topography_surface)
    set_topography_boundary!(model, topography_surface)
    # Compute the stretched vertical grid.
    for kz in (-nbz):(nz+nbz)
        ztildes[kz] = map(z[kz] + 0.5 * dz, lz)
    end
    for kz in (-nbz+1):(nz+nbz)
        zs[kz] = 0.5 * (ztildes[kz] + ztildes[kz-1])
    end
    zs[-nbz] = ztildes[-nbz] - 0.5 * (ztildes[nbz+1] - ztildes[nbz])

    # Compute the physical layers.
    for kz in (-nbz):(nz+nbz)
        ztfc[:, :, kz] .= (lz[1] .- topography_surface) / lz[1] * zs[kz] .+
                          topography_surface
    end
    end # timer
end
function initialize_atmosphere!(model::Model)
    @trixi_timeit timer() "Initialize atmosphere" begin
    #! format: noindent
    (; nx, ny, nz) = model.domain
    (; nbx, nby, nbz) = model.parameters.domain
    (; jac, ztfc, dz) = model.grid
    # (; pStrat, rhoStrat, thetaStrat, bvsStrat, jac) = cache
    (; pstrattfc, rhostrattfc, thetastrattfc, bvsstrattfc, t0, p0) = model.atmosphere
    # (; gamma, gamma_1, kappa, kappaInv, gammaInv, Rsp, g, rhoRef, pRef, aRef, uRef, lRef, tRef, thetaRef, Ma, Fr, kappa, sig, press0_dim, Temp0_dim, T0, N2, NN, mu_viscous_dim, ReInv, Re,) = equations

    setup_topography!(model)

    c = model.constants
    for ix in (-nbx):(nx+nbx), jy in (-nby):(ny+nby)
        for kz in -1:(nz+2)
            pstrattfc[ix, jy, kz] = p0 * exp(-c.sig * ztfc[ix, jy, kz] / c.gamma / t0)
            thetastrattfc[ix, jy, kz] = t0 * exp(c.kappa * c.sig / t0 * ztfc[ix, jy, kz])
            rhostrattfc[ix, jy, kz] = pstrattfc[ix, jy, kz] / thetastrattfc[ix, jy, kz]
        end
    end

    for ix in (-nbx):(nx+nbx), jy in (-nby):(ny+nby)
        bvsstrattfc[ix, jy, -1] = c.g_ndim / thetastrattfc[ix, jy, 0] / jac(ix, jy, 0) *
                                  (thetastrattfc[ix, jy, 1] - thetastrattfc[ix, jy, 0]) / dz
        for kz in 1:nz
            bvsstrattfc[ix, jy, kz] = c.g_ndim / thetastrattfc[ix, jy, kz] / jac(ix, jy, kz) *
                                      0.5 *
                                      (thetastrattfc[ix, jy, kz+1] -
                                       thetastrattfc[ix, jy, kz-1]) /
                                      dz
        end
        bvsstrattfc[ix, jy, nz+1] = c.g_ndim / thetastrattfc[ix, jy, nz+1] /
                                    jac(ix, jy, nz + 1) *
                                    (thetastrattfc[ix, jy, nz+1] - thetastrattfc[ix, jy, nz]) /
                                    dz
        bvsstrattfc[ix, jy, nz+2] = bvsstrattfc[ix, jy, nz+1]
    end
    end # timer
end
