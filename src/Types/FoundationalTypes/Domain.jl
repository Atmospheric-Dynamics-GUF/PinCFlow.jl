struct Domain{
    A <: MPI.Comm,
    B <: Bool,
    C <: Integer,
    D <: AbstractVector{<:AbstractFloat},
}

    # MPI variables.
    comm::A
    master::B
    rank::C
    root::C

    # Local grid size.
    nx::C
    ny::C
    nz::C

    # Local grid size with boundary cells
    nxx::C
    nyy::C
    nzz::C

    # Global grid size with boundary cells.
    sizexx::C
    sizeyy::C
    sizezz::C

    # Index offsets.
    io::C
    jo::C
    ko::C

    # Index bounds.
    i0::C
    i1::C
    j0::C
    j1::C
    k0::C
    k1::C

    # Source and destination ranks for halos.
    left::C
    right::C
    backward::C
    forward::C
    down::C
    up::C

    # Auxiliary arrays for horizontal averages.
    local_sum::D
    global_sum::D
end

function Domain(namelists::Namelists)
    (; sizex, sizey, sizez, nbx, nby, nbz, npx, npy, npz) = namelists.domain
    (; zboundaries) = namelists.setting

    # Initialize MPI.
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    root = 0
    if rank == root
        master = true
    else
        master = false
    end
    np = MPI.Comm_size(MPI.COMM_WORLD)

    # Check if parallelization is set up correctly.
    if master && npx * npy * npz != np
        error("Error in Domain: npx * npy * npz != np!")
    end
    if master && npx > 1 && nbx > div(sizex, npx)
        error("Error in Domain: npx > 1 && nbx > div(sizex, npx)!")
    end
    if master && npy > 1 && nby > div(sizey, npy)
        error("Error in Domain: npy > 1 && nby > div(sizey, npy)!")
    end
    if master && npz > 1 && nbz > div(sizez, npz)
        error("Error in Domain: npz > 1 && nbz > div(sizez, npz)!")
    end

    # Set dimensions and periodicity.
    dims = [npx, npy, npz]
    if zboundaries == SolidWallBoundaries()
        periods = [true, true, false]
    else
        error("Error in Domain: Unknown zboundaries!")
    end

    # Create a Cartesian topology.
    comm = MPI.Cart_create(MPI.COMM_WORLD, dims; periodic = periods)
    rank = MPI.Comm_rank(comm)
    coords = MPI.Cart_coords(comm, rank)

    # Set local grid size.
    if coords[1] == npx - 1
        nx = div(sizex, npx) + sizex % npx
    else
        nx = div(sizex, npx)
    end
    if coords[2] == npy - 1
        ny = div(sizey, npy) + sizey % npy
    else
        ny = div(sizey, npy)
    end
    if coords[3] == npz - 1
        nz = div(sizez, npz) + sizez % npz
    else
        nz = div(sizez, npz)
    end

    # Set grid sizes with boundary cells.
    nxx = nx + 2 * nbx
    nyy = ny + 2 * nby
    nzz = nz + 2 * nbz
    sizexx = sizex + 2 * nbx
    sizeyy = sizey + 2 * nby
    sizezz = sizez + 2 * nbz

    # Set index offsets.
    io = coords[1] * div(sizex, npx)
    jo = coords[2] * div(sizey, npy)
    ko = coords[3] * div(sizez, npz)

    # Set index bounds.
    i0 = nbx + 1
    i1 = i0 + nx - 1
    j0 = nby + 1
    j1 = j0 + ny - 1
    k0 = nbz + 1
    k1 = k0 + nz - 1

    # Find the neighbour processors.
    (left, right) = MPI.Cart_shift(comm, 0, 1)
    (backward, forward) = MPI.Cart_shift(comm, 1, 1)
    (down, up) = MPI.Cart_shift(comm, 2, 1)

    # Initialize auxiliary arrays for horizontal averages.
    (local_sum, global_sum) = (zeros(sizez) for i in 1:2)

    # Return Domain instance.
    return Domain(
        comm,
        master,
        rank,
        root,
        nx,
        ny,
        nz,
        nxx,
        nyy,
        nzz,
        sizexx,
        sizeyy,
        sizezz,
        io,
        jo,
        ko,
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
        left,
        right,
        backward,
        forward,
        down,
        up,
        local_sum,
        global_sum,
    )
end
