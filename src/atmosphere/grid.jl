using OffsetArrays

struct Jacobian
    topography_surface::Any
    zTildeS::Any
    dz::Any
    lz::Any
end

struct MetricTensor
    topography_surface::Any
    zTildeS::Any
    dx::Any
    dy::Any
    dz::Any
    lz::Any
    zS::Any
end

struct Grid{F <: AbstractFloat,
            VOF <: OffsetVector{F},
            MOF <: OffsetMatrix{F},
            A3OF <: OffsetArray{F, 3}}
    # A5OF<:OffsetArray{F,5}}
    # Scaled domain.
    lx::VOF
    ly::VOF
    lz::VOF
    # Grid spacings.
    dx::F
    dy::F
    dz::F
    # Coordinates.
    x::VOF
    y::VOF
    z::VOF
    # Stretched vertical grid.
    zs::VOF
    ztildes::VOF
    # Topography.
    topography_surface::MOF
    # Jacobian and metric tensor.
    jac::Jacobian
    met::MetricTensor
    # Vertical layers.
    ztfc::A3OF
    ztildetfc::A3OF
    # Boundaries
    xboundary::Any
    yboundary::Any
    zboundary::Any
end

function Grid(p::Parameters, lRef) # todo: get rid of lref
    nx, ny, nz = p.domain.sizex, p.domain.sizey, p.domain.sizez
    nbx, nby, nbz = p.domain.nbx, p.domain.nby, p.domain.nbz
    topography_surface = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby),
                                     (-nbx):(nx + nbx), (-nby):(ny + nby))
    zTildeTFC = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby, nz + 1 + 2 * nbz),
                            (-nbx):(nx + nbx),
                            (-nby):(ny + nby),
                            (-nbz):(nz + nbz))
    zTFC = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby, nz + 1 + 2 * nbz),
                       (-nbx):(nx + nbx),
                       (-nby):(ny + nby),
                       (-nbz):(nz + nbz))
    zTildeS = OffsetArray(zeros(nz + 1 + 2 * nbz), (-nbz):(nz + nbz))
    zS = OffsetArray(zeros(nz + 1 + 2 * nbz), (-nbz):(nz + nbz))

    lx_dim = OffsetArray(p.domain.lx_dim, 0:1)
    ly_dim = OffsetArray(p.domain.ly_dim, 0:1)
    lz_dim = OffsetArray(p.domain.lz_dim, 0:1)

    lx = lx_dim ./ lRef # TODO - Keep lx_dim in `grid` and lx_ref in semi.physics
    ly = ly_dim ./ lRef
    lz = lz_dim ./ lRef

    dx = (lx[1] - lx[0]) / nx
    dy = (ly[1] - ly[0]) / ny
    dz = (lz[1] - lz[0]) / nz

    sizeX = nx
    sizeY = ny
    sizeZ = nz

    x = OffsetArray(zeros(sizeX + 1 + 2nbx), (-nbx):(sizeX + nbx))
    y = OffsetArray(zeros(sizeY + 1 + 2nby), (-nby):(sizeY + nby))
    z = OffsetArray(zeros(sizeZ + 1 + 2nbz), (-nbz):(sizeZ + nbz))

    for i in (-nbx):(sizeX + nbx)
        x[i] = lx[0] + (i - 1) * dx + dx / 2.0
    end

    for j in (-nby):(sizeY + nby)
        y[j] = ly[0] + (j - 1) * dy + dy / 2.0
    end

    for k in (-nbz):(sizeZ + nbz)
        z[k] = lz[0] + (k - 1) * dz + dz / 2.0
    end

    # Build arrays.
    jac = Jacobian(topography_surface, zTildeS, dz, lz)
    met = MetricTensor(topography_surface, zTildeS, dx, dy, dz, lz, zS)

    xboundary = BoundaryCondition(p.boundaries.xboundary)
    yboundary = BoundaryCondition(p.boundaries.yboundary)
    zboundary = BoundaryCondition(p.boundaries.zboundary)
    # Return constructed Grid
    return Grid(lx, ly, lz, dx, dy, dz, x, y, z, zS, zTildeS, topography_surface, jac, met,
                zTFC, zTildeTFC, xboundary, yboundary, zboundary)
end

function BoundaryCondition(s::AbstractString)
    if s == "periodic"
        return PeriodicBC()
    elseif s == "solid_wall"
        return SolidWallBC()
    else
        error("unsupported boundary condition $s")
    end
end

function (met::MetricTensor)(i, j, k, mu, nu)
    # Metric tensor.
    (; topography_surface, zS, zTildeS, dx, dy, dz, lz) = met
    if (mu == 1 && nu == 3) || (mu == 3 && nu == 1)
        return (topography_surface[i + 1, j] - topography_surface[i - 1, j]) / (2.0 * dx) *
               (zS[k] - lz[1]) / (lz[1] - topography_surface[i, j]) * dz /
               (zTildeS[k] - zTildeS[k - 1])
    elseif (mu == 2 && nu == 3) || (mu == 3 && nu == 2)
        return (topography_surface[i, j + 1] - topography_surface[i, j - 1]) / (2.0 * dy) *
               (zS[k] - lz[1]) / (lz[1] - topography_surface[i, j]) * dz /
               (zTildeS[k] - zTildeS[k - 1])
    elseif mu == 3 && nu == 3
        return ((lz[1] / (lz[1] - topography_surface[i, j]))^2.0 +
                ((zS[k] - lz[1]) / (lz[1] - topography_surface[i, j]))^2.0 *
                (((topography_surface[i + 1, j] - topography_surface[i - 1, j]) /
                  (2.0 * dx))^2.0 +
                 ((topography_surface[i, j + 1] - topography_surface[i, j - 1]) /
                  (2.0 * dy))^2.0)) * (dz / (zTildeS[k] - zTildeS[k - 1]))^2.0
    else
        @assert false "UNDEFINED CASE!!!"
    end
end

function (jac::Jacobian)(i, j, k)
    # Jacobian.
    (; topography_surface, zTildeS, dz, lz) = jac
    return (lz[1] - topography_surface[i, j]) / lz[1] * (zTildeS[k] - zTildeS[k - 1]) / dz
end

function Base.getindex(jac::Jacobian, i, j, k)
    return jac(i, j, k)
end

function Base.getindex(met::MetricTensor, i, j, k, u, v)
    return met(i, j, k, u, v)
end
