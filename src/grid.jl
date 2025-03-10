using OffsetArrays

struct Grid{F<:AbstractFloat,
    VOF<:OffsetVector{F},
    MOF<:OffsetMatrix{F},
    A3OF<:OffsetArray{F,3},
    A5OF<:OffsetArray{F,5}}
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
    jac::A3OF
    met::A5OF
    # Vertical layers.
    ztfc::A3OF
    ztildetfc::A3OF
    # Boundaries
    xboundary
    yboundary
    zboundary
end

function Grid(p::Parameters, lRef) # todo: get rid of lref
    nx, ny, nz = p.domain.sizex, p.domain.sizey, p.domain.sizez
    nbx, nby, nbz = p.domain.nbx, p.domain.nby, p.domain.nbz
    topography_surface = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby),
        (-nbx):(nx+nbx), (-nby):(ny+nby))
    zTildeTFC = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby, nz + 1 + 2 * nbz),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz))
    zTFC = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby, nz + 1 + 2 * nbz),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz))
    zTildeS = OffsetArray(zeros(nz + 1 + 2 * nbz), (-nbz):(nz+nbz))
    zS = OffsetArray(zeros(nz + 1 + 2 * nbz), (-nbz):(nz+nbz))

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

    x = OffsetArray(zeros(sizeX + 1 + 2nbx), (-nbx):(sizeX+nbx))
    y = OffsetArray(zeros(sizeY + 1 + 2nby), (-nby):(sizeY+nby))
    z = OffsetArray(zeros(sizeZ + 1 + 2nbz), (-nbz):(sizeZ+nbz))

    for i in (-nbx):(sizeX+nbx)
        x[i] = lx[0] + (i - 1) * dx + dx / 2.0
    end

    for j in (-nby):(sizeY+nby)
        y[j] = ly[0] + (j - 1) * dy + dy / 2.0
    end

    for k in (-nbz):(sizeZ+nbz)
        z[k] = lz[0] + (k - 1) * dz + dz / 2.0
    end
    dx = p.domain.sizex / p.domain.nbx
    dy = p.domain.sizey / p.domain.nby
    dz = p.domain.sizez / p.domain.nbz


    # Build arrays.
    jac = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby, nz + 1 + 2 * nbz),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz))

    met = OffsetArray(zeros(nx + 1 + 2 * nbx, ny + 1 + 2 * nby, nz + 1 + 2 * nbz, 3, 3),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz),
        1:3,
        1:3)

    xboundary = PeriodicBC()
    yboundary = PeriodicBC()
    zboundary = PeriodicBC()
    # Return constructed Grid
    return Grid(lx, ly, lz, dx, dy, dz, x, y, z, zS, zTildeS, topography_surface, jac, met, zTFC, zTildeTFC, xboundary, yboundary, zboundary)
end
