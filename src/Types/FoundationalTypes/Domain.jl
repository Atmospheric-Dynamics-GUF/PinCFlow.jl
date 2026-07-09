"""
```julia
Domain{A <: MPI.Comm, B <: Bool, C <: Integer}
```

Collection of domain-decomposition and MPI-communication parameters.

```julia
Domain(namelists::Namelists; base_comm::MPI.Comm = MPI.COMM_WORLD)::Domain
```

Construct a `Domain` instance from the model parameters in `namelists`.

This method first creates a Cartesian topology from the base communicator `base_comm`, with periodic boundaries in the first two dimensions (``\\hat{x}`` and ``\\hat{y}``) but not in the last (``\\hat{z}``). The domain is divided into corresponding subdomains, where in each direction, the number of grid points (`nx`, `ny` and `nz`) is the result of floor division of the global grid size (`namelists.domain.x_size`, `namelists.domain.y_size` and `namelists.domain.z_size`) by the number of processes in that direction (`namelists.domain.npx`, `namelists.domain.npy` and `namelists.domain.npz`). The remainder of the floor division is included in the grid-point count of the last processes (in each direction). The index bounds (`(i0, i1)`, `(j0, j1)` and `(k0, k1)`) are set such that they exclude the first and last `namelists.domain.nbx`, `namelists.domain.nby` and `namelists.domain.nbz` cells in ``\\hat{x}``, ``\\hat{y}`` and ``\\hat{z}``, respectively (these are not included in `nx`, `ny` and `nz`).

# Fields

General MPI communication:

  - `comm::A`: MPI communicator with Cartesian topology for the computational domain.

  - `master::B`: Boolean flag indicating if this process is the master process (rank 0).

  - `rank::C`: MPI rank of this process within the communicator `comm`.

  - `root::C`: Root process rank (0).

Dimensions of the MPI subdomain:

  - `nx::C`: Number of physical grid points in ``\\hat{x}``-direction.

  - `ny::C`: Number of physical grid points in ``\\hat{y}``-direction.

  - `nz::C`: Number of physical grid points in ``\\hat{z}``-direction.

  - `nxx::C`: Number of computational grid points in ``\\hat{x}``-direction (including halo/boundary cells).

  - `nyy::C`: Number of computational grid points in ``\\hat{y}``-direction (including halo/boundary cells).

  - `nzz::C`: Number of computational grid points in ``\\hat{z}``-direction (including halo/boundary cells).

Index offsets and bounds:

  - `io::C`: MPI offset in ``\\hat{x}``-direction.

  - `jo::C`: MPI offset in ``\\hat{y}``-direction.

  - `ko::C`: MPI offset in ``\\hat{z}``-direction.

  - `i0::C`: First physical grid cell of the subdomain in ``\\hat{x}``-direction.

  - `i1::C`: Last physical grid cell of the subdomain in ``\\hat{x}``-direction.

  - `j0::C`: First physical grid cell of the subdomain in ``\\hat{y}``-direction.

  - `j1::C`: Last physical grid cell of the subdomain in ``\\hat{y}``-direction.

  - `k0::C`: First physical grid cell of the subdomain in ``\\hat{z}``-direction.

  - `k1::C`: Last physical grid cell of the subdomain in ``\\hat{z}``-direction.

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

# Keywords

  - `base_comm`: MPI communicator which is used to create the Cartesian communicator for the integration.
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

@ivy function Domain(
    namelists::Namelists;
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::Domain
    (;
        x_size,
        y_size,
        z_size,
        nbx,
        nby,
        nbz,
        npx,
        npy,
        npz,
        vertical_boundary_condition,
    ) = namelists.domain

    # Initialize MPI.
    !MPI.Initialized() && MPI.Init()
    rank = MPI.Comm_rank(base_comm)
    root = 0
    if rank == root
        master = true
    else
        master = false
    end
    np = MPI.Comm_size(base_comm)

    # Check if there will be enough boundary cells.
    if master && nbx < 3
        error("Too few boundary cells: nbx < 3!")
    end
    if master && nby < 3
        error("Too few boundary cells: nby < 3!")
    end
    if master && nbz < 3
        error("Too few boundary cells: nbz < 3!")
    end

    # Check if parallelization is set up correctly.
    if master && npx * npy * npz != np
        error("Incorrect domain decomposition: npx * npy * npz != np!")
    end
    if master && x_size > 1 && nbx > div(x_size, npx)
        error("Too many MPI subdomains: x_size > 1 && nbx > div(x_size, npx)!")
    end
    if master && y_size > 1 && nby > div(y_size, npy)
        error("Too many MPI subdomains: y_size > 1 && nby > div(y_size, npy)!")
    end
    if master && z_size > 1 && nbz > div(z_size, npz)
        error("Too many MPI subdomains: z_size > 1 && nbz > div(z_size, npz)!")
    end

    # Set dimensions and periodicity.
    dims = [npx, npy, npz]

    if vertical_boundary_condition == :SolidWall
        periods = [true, true, false]
    elseif vertical_boundary_condition == :Periodic
        periods = [true, true, true]
    else 
        error("Incorrect choice of `vertical_boundary_condition`. Must be `:SolidWall` or `:Periodic`")
    end

    # Create a Cartesian topology.
    comm = MPI.Cart_create(base_comm, dims; periodic = periods)
    rank = MPI.Comm_rank(comm)
    coords = MPI.Cart_coords(comm, rank)

    # Set local grid size.
    if coords[1] == npx - 1
        nx = div(x_size, npx) + x_size % npx
    else
        nx = div(x_size, npx)
    end
    if coords[2] == npy - 1
        ny = div(y_size, npy) + y_size % npy
    else
        ny = div(y_size, npy)
    end
    if coords[3] == npz - 1
        nz = div(z_size, npz) + z_size % npz
    else
        nz = div(z_size, npz)
    end

    # Set grid sizes with boundary cells.
    nxx = nx + 2 * nbx
    nyy = ny + 2 * nby
    nzz = nz + 2 * nbz

    # Set index offsets.
    io = coords[1] * div(x_size, npx)
    jo = coords[2] * div(y_size, npy)
    ko = coords[3] * div(z_size, npz)

    # Set index bounds.
    i0 = nbx + 1
    i1 = i0 + nx - 1
    j0 = nby + 1
    j1 = j0 + ny - 1
    k0 = nbz + 1
    k1 = k0 + nz - 1

    # Find the neighbor processors.
    (left, right) = MPI.Cart_shift(comm, 0, 1)
    (backward, forward) = MPI.Cart_shift(comm, 1, 1)
    (down, up) = MPI.Cart_shift(comm, 2, 1)

    # Create communicators for horizontal and vertical averages.
    layer_comm = MPI.Comm_split(comm, coords[3], rank)
    column_comm = MPI.Comm_split(comm, coords[2] * npx + coords[1], rank)

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
