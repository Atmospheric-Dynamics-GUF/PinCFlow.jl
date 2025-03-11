
"""
Properties of the MPI domain.
"""
struct Domain{C,B<:Bool,I<:Integer}

    # MPI variables.
    comm::C
    master::B
    rank::I
    root::I

    # Number of CPUs.
    idim::I
    jdim::I

    # Start indices of the local grid.
    is::I
    js::I

    # Local grid size.
    nx::I
    ny::I
    nz::I

    # halo cells
    nbx::I
    nby::I
    nbz::I
    # Local grid size including halo/ghost cells.
    nxx::I
    nyy::I
    nzz::I
end

function init_mpi()
    (nothing, false, 1, 1)
end
function Domain(params::Parameters)
    # TODO
    comm, master, rank, root = init_mpi()
    idim = params.domain.nprocx
    jdim = params.domain.nprocy
    (; nbx, nby, nbz) = params.domain
    js, is = 1, 1
    nx, ny, nz = params.domain.sizex, params.domain.sizey, params.domain.sizez
    nxx = nx + 2 * nbx + 1 # TODO: +1 or not?
    nyy = ny + 2 * nby + 1
    nzz = nz + 2 * nbz + 1

    return Domain(comm, master, rank, root, idim, jdim, js, is,
        nx, ny, nz, nbx, nby, nbz, nxx, nyy, nzz)
end
