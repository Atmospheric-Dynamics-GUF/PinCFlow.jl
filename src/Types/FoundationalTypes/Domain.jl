"""
```julia
Domain{A <: MPI.Comm, B <: Bool, C <: Integer}
```

Collection of domain-decomposition and MPI-communication parameters.

```julia
Domain(namelists::Namelists)
```

Construct a `Domain` instance from the model parameters in `namelists`.

If `namelists.domain.base_comm` is equal to `MPI.COMM_WORLD`, this method first initializes the MPI parallelization by calling `MPI.Init()`. It then creates a Cartesian topology from the base communicator, with periodic boundaries in the first two dimensions (``\\widehat{x}`` and ``\\widehat{y}``) but not in the last (``\\widehat{z}``). The domain is divided into corresponding subdomains, where in each direction, the number of grid points (`nx`, `ny` and `nz`) is the result of floor division of the global grid size (`namelists.domain.sizex`, `namelists.domain.sizey` and `namelists.domain.sizez`) by the number of processes in that direction (`npx`, `npy` and `npz`). The remainder of the floor division is included in the grid-point count of the last processes (in each direction). The index bounds (`(i0, i1)`, `(j0, j1)` and `(k0, k1)`) are set such that they exclude the first and last `namelists.domain.nbx`, `namelists.domain.nby` and `namelists.domain.nbz` cells in ``\\widehat{x}``, ``\\widehat{y}`` and ``\\widehat{z}``, respectively.

# Fields

General MPI communication:

- `comm::A`: MPI communicator with Cartesian topology for the computational domain.
- `master::B`: Boolean flag indicating if this process is the master process (rank 0).
- `rank::C`: MPI rank of this process within the communicator `comm`.
- `root::C`: Root process rank (0).

Dimensions of the MPI subdomain:

- `nx::C`: Number of physical grid points in ``\\widehat{x}``-direction.
- `ny::C`: Number of physical grid points in ``\\widehat{y}``-direction.
- `nz::C`: Number of physical grid points in ``\\widehat{z}``-direction.
- `nx::C`: Number of computational grid points in ``\\widehat{x}``-direction (including halo/boundary cells).
- `ny::C`: Number of computational grid points in ``\\widehat{y}``-direction (including halo/boundary cells).
- `nz::C`: Number of computational grid points in ``\\widehat{z}``-direction (including halo/boundary cells).

Dimensions of the entire domain:

- `sizexx::C`: Number of computational grid points in ``\\widehat{x}``-direction (including halo/boundary cells).
- `sizeyy::C`: Number of computational grid points in ``\\widehat{y}``-direction (including halo/boundary cells).
- `sizezz::C`: Number of computational grid points in ``\\widehat{z}``-direction (including halo/boundary cells).

Index offsets and bounds:

- `io::C`: MPI offset in ``\\widehat{x}``-direction.
- `jo::C`: MPI offset in ``\\widehat{y}``-direction.
- `ko::C`: MPI offset in ``\\widehat{z}``-direction.
- `i0::C`: First physical grid cell of the subdomain in ``\\widehat{x}``-direction.
- `i1::C`: Last physical grid cell of the subdomain in ``\\widehat{x}``-direction.
- `j0::C`: First physical grid cell of the subdomain in ``\\widehat{y}``-direction.
- `j1::C`: Last physical grid cell of the subdomain in ``\\widehat{y}``-direction.
- `k0::C`: First physical grid cell of the subdomain in ``\\widehat{z}``-direction.
- `k1::C`: Last physical grid cell of the subdomain in ``\\widehat{z}``-direction.

Neighbor-process ranks:

- `left::C`: Rank of the next process to the left (negative ``x``-direction).
- `right::C`: Rank of the next process to the right (positive ``x``-direction).
- `backward::C`: Rank of the next process to the back (negative ``y``-direction).
- `forward::C`: Rank of the next process to the front (positive ``y``-direction).
- `down::C`: Rank of the next process to the bottom (negative ``z``-direction).
- `up::C`: Rank of the next process to the top (positive ``z``-direction).

Horizontal and vertical communication:

- `layer_comm::A`: MPI communicator for processes in the same layer (i.e. with the same vertical index).
- `column_comm::A`: MPI communicator for processes in the same column (i.e. with the same horizontal indices).

# Arguments

- `namelists`: Namelists with all model parameters.
"""
struct Domain{A <: MPI.Comm, B <: Bool, C <: Integer}

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

    # Communicators for horizontal and vertical averages.
    layer_comm::A
    column_comm::A
end

function Domain(namelists::Namelists)
    (; sizex, sizey, sizez, nbx, nby, nbz, npx, npy, npz, base_comm) =
        namelists.domain
    (; zboundaries) = namelists.setting

    # Initialize MPI.
    if base_comm == MPI.COMM_WORLD
        MPI.Init()
    end
    rank = MPI.Comm_rank(base_comm)
    root = 0
    if rank == root
        master = true
    else
        master = false
    end
    np = MPI.Comm_size(base_comm)

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
    comm = MPI.Cart_create(base_comm, dims; periodic = periods)
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

    # Create communicators for horizontal and vertical averages.
    layer_comm = MPI.Comm_split(comm, coords[3], rank)
    column_comm = MPI.Comm_split(comm, coords[2] * npx + coords[1], rank)

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
        layer_comm,
        column_comm,
    )
end
