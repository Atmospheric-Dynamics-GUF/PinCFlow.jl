"""
```julia
Domain{A <: MPI.Comm, B <: Bool, C <: Integer}
```

Collection of domain-decomposition and MPI-communication parameters.

# Fields

## MPI Communication

  - `comm`: MPI communicator with Cartesian topology for the computational domain
  - `master`: Boolean flag indicating if this process is the master process (rank 0)
  - `rank`: MPI rank of this process within the communicator
  - `root`: Root process rank (typically 0)

## Local Grid Dimensions

  - `nx`, `ny`, `nz`: Number of grid points in each direction for this process (excluding boundary cells)
  - `nxx`, `nyy`, `nzz`: Total grid points in each direction including boundary cells (`n* + 2*nb*`)

## Global Grid Dimensions

  - `sizexx`, `sizeyy`, `sizezz`: Global grid size in each direction including boundary cells

## Index Management

  - `io`, `jo`, `ko`: Index offsets in each direction for global-to-local coordinate transformation
  - `i0`, `i1`: Starting and ending indices for computational domain in x-direction
  - `j0`, `j1`: Starting and ending indices for computational domain in y-direction
  - `k0`, `k1`: Starting and ending indices for computational domain in z-direction

## Neighbor Process Ranks

  - `left`, `right`: MPI ranks of neighbor processes in x-direction (negative/positive)
  - `backward`, `forward`: MPI ranks of neighbor processes in y-direction (negative/positive)
  - `down`, `up`: MPI ranks of neighbor processes in z-direction (negative/positive)

## Specialized Communicators

  - `layer_comm`: MPI communicator for processes in the same horizontal layer (same k-coordinate)
  - `column_comm`: MPI communicator for processes in the same vertical column (same i,j coordinates)

# Usage

Arrays should be allocated with dimensions `[nxx, nyy, nzz]` and accessed using
indices `i0:i1`, `j0:j1`, `k0:k1` for the computational domain.
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

"""
```julia
Domain(namelists::Namelists) -> Domain
```

Construct a `Domain` instance for MPI domain decomposition.

This constructor initializes the MPI environment and sets up a 3D Cartesian topology
for parallel computation. It distributes the global computational domain across
multiple MPI processes and establishes communication patterns for halo exchanges.

# Arguments

  - `namelists`: Configuration object containing domain and setting parameters

# Returns

  - `::Domain`: `Domain` instance.

# Error Conditions

Throws an error if:

  - Process decomposition doesn't match total process count (`npx * npy * npz != np`)
  - Boundary cell count exceeds local domain size in any parallelized direction
  - Unknown boundary condition type is specified

Currently supports `SolidWallBoundaries()` for z-direction with periodic boundaries in x and y.
"""
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
