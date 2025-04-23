struct Domain{
    A <: MPI.Comm,
    B <: Bool,
    C <: Integer,
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 5},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 5},
    H <: AbstractMatrix{<:AbstractFloat},
    I <: AbstractMatrix{<:AbstractFloat},
    J <: AbstractVector{<:AbstractFloat},
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
    back::C
    forw::C

    # Auxiliary arrays for setting halos of a field.
    send_f3_left::D
    send_f3_right::D
    recv_f3_left::D
    recv_f3_right::D
    send_f5_left::E
    send_f5_right::E
    recv_f5_left::E
    recv_f5_right::E
    send_f3_back::F
    send_f3_forw::F
    recv_f3_back::F
    recv_f3_forw::F
    send_f5_back::G
    send_f5_forw::G
    recv_f5_back::G
    recv_f5_forw::G

    # Auxiliary arrays for setting halos of a reduced field.
    send_rf3_left::H
    send_rf3_right::H
    recv_rf3_left::H
    recv_rf3_right::H
    send_rf3_back::I
    send_rf3_forw::I
    recv_rf3_back::I
    recv_rf3_forw::I

    # Auxiliary arrays for horizontal averaging.
    local_sum::J
    global_sum::J
end

function Domain(namelists::Namelists)

    # Get domain parameters.
    (; sizex, sizey, sizez, nbx, nby, nbz, nprocx, nprocy) = namelists.domain

    # Initialize MPI.
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    root = 0
    if rank == root
        master = true
    else
        master = false
    end
    nproc = MPI.Comm_size(MPI.COMM_WORLD)

    # Check if parallelization is set up correctly.
    if master && nprocx * nprocy != nproc
        error("Error in Domain: nprocx * nprocy != nproc!")
    end
    if master && nprocx > 1 && nbx > div(sizex, nprocx)
        error("Error in Domain: nprocx > 1 && nbx > div(sizex, nprocx)!")
    end
    if master && nprocy > 1 && nby > div(sizey, nprocy)
        error("Error in Domain: nprocy > 1 && nby > div(sizey, nprocy)!")
    end

    # Set dimensions and periodicity.
    dims = [nprocx, nprocy]
    periods = [true, true]

    # Create a Cartesian topology.
    comm = MPI.Cart_create(MPI.COMM_WORLD, dims; periodic = periods)
    rank = MPI.Comm_rank(comm)
    coords = MPI.Cart_coords(comm, rank)

    # Set local grid size.
    if coords[1] == nprocx - 1
        nx = div(sizex, nprocx) + sizex % nprocx
    else
        nx = div(sizex, nprocx)
    end
    if coords[2] == nprocy - 1
        ny = div(sizey, nprocy) + sizey % nprocy
    else
        ny = div(sizey, nprocy)
    end
    nz = sizez

    # Set grid sizes with boundary cells.
    nxx = nx + 2 * nbx
    nyy = ny + 2 * nby
    nzz = nz + 2 * nbz
    sizexx = sizex + 2 * nbx
    sizeyy = sizey + 2 * nby
    sizezz = sizez + 2 * nbz

    # Set index offsets.
    io = coords[1] * div(sizex, nprocx)
    jo = coords[2] * div(sizey, nprocy)

    # Set index bounds.
    i0 = nbx + 1
    i1 = i0 + nx - 1
    j0 = nby + 1
    j1 = j0 + ny - 1
    k0 = nbz + 1
    k1 = k0 + nz - 1

    # Find the neighbour processors.
    (left, right) = MPI.Cart_shift(comm, 0, 1)
    (back, forw) = MPI.Cart_shift(comm, 1, 1)

    # Initialize auxiliary arrays for setting all halo layers.
    (send_f3_left, send_f3_right, recv_f3_left, recv_f3_right) =
        (zeros(nbx, nyy, nzz) for i in 1:4)
    (send_f5_left, send_f5_right, recv_f5_left, recv_f5_right) =
        (zeros(nbx, nyy, nzz, 3, 2) for i in 1:4)
    (send_f3_back, send_f3_forw, recv_f3_back, recv_f3_forw) =
        (zeros(nxx, nby, nzz) for i in 1:4)
    (send_f5_back, send_f5_forw, recv_f5_back, recv_f5_forw) =
        (zeros(nxx, nby, nzz, 3, 2) for i in 1:4)

    # Initialize auxiliary arrays for setting one halo layer.
    (send_rf3_left, send_rf3_right, recv_rf3_left, recv_rf3_right) =
        (zeros(ny + 2, nz + 2) for i in 1:4)
    (send_rf3_back, send_rf3_forw, recv_rf3_back, recv_rf3_forw) =
        (zeros(nx + 2, nz + 2) for i in 1:4)

    # Initialize auxiliary arrays for horizontal averaging.
    (local_sum, global_sum) = (zeros(nz) for i in 1:2)

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
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
        left,
        right,
        back,
        forw,
        send_f3_left,
        send_f3_right,
        recv_f3_left,
        recv_f3_right,
        send_f5_left,
        send_f5_right,
        recv_f5_left,
        recv_f5_right,
        send_f3_back,
        send_f3_forw,
        recv_f3_back,
        recv_f3_forw,
        send_f5_back,
        send_f5_forw,
        recv_f5_back,
        recv_f5_forw,
        send_rf3_left,
        send_rf3_right,
        recv_rf3_left,
        recv_rf3_right,
        send_rf3_back,
        send_rf3_forw,
        recv_rf3_back,
        recv_rf3_forw,
        local_sum,
        global_sum,
    )
end
