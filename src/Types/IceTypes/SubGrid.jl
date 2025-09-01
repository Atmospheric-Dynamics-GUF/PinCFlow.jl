struct SubGrid{A <: AbstractFloat, B <: Integer} 
    dxsc :: A
    dysc :: A
    dzsc :: A

    nxnscx :: B
    nynscy :: B
    nznscz :: B

    nxnscxx :: B
    nynscyy :: B
    nznsczz :: B

end

function SubGrid(namelists::Namelists, domain :: Domain, grid::Grid)

    (; icesetup) = namelists.ice
    return SubGrid(namelists, domain, grid, icesetup)
end

function SubGrid(namelists::Namelists, 
    domain :: Domain, 
    grid :: Grid,
    icesetup::NoIce
    )

    (; nbx, nby, nbz) = namelists.domain
    (; nx, ny, nz) = domain
    (; dx, dy, dz) = grid

    dxsc = dx
    dysc = dy
    dzsc = dz

    nxnscx = nx
    nynscy = ny
    nznscz = nz

    nxnscxx = nxnscx + 2 * nbx
    nynscyy = nynscy + 2 * nby
    nznsczz = nznscz + 2 * nbz

    return SubGrid(dxsc, dysc, dzsc, 
    nxnscx, nynscy, nznscz, 
    nxnscxx, nynscyy, nznsczz)

end    

function SubGrid(namelists::Namelists, 
    domain :: Domain,
    grid::Grid,
    icesetup::AbstractIce
    )

    (; compute_cloudcover, nscx, nscy, nscz) = namelists.ice
    (; nbx, nby, nbz) = namelists.domain
    (; nx, ny, nz) = domain
    (; dx, dy, dz) = grid

    if compute_cloudcover == 2

        dxsc = dx / nscx
        dysc = dy / nscy
        dzsc = dz / nscz

        nxnscx = nx * nscx
        nynscy = ny * nscy
        nznscz = nz * nscz

        nxnscxx = nxnscx + 2 * nbx * nscx
        nynscyy = nynscy + 2 * nby * nscy
        nznsczz = nznscz + 2 * nbz * nscz

    else
        dxsc = dx
        dysc = dy
        dzsc = dz

        nxnscx = nx
        nynscy = ny
        nznscz = nz

        nxnscxx = nxnscx + 2 * nbx
        nynscyy = nynscy + 2 * nby
        nznsczz = nznscz + 2 * nbz
    end

    return SubGrid(dxsc, dysc, dzsc, 
    nxnscx, nynscy, nznscz, 
    nxnscxx, nynscyy, nznsczz)

end