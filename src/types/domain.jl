struct Domain{
  A <: MPI.Comm,
  B <: Bool,
  C <: Integer,
  D <: AbstractMatrix{<:AbstractFloat},
  E <: AbstractArray{<:AbstractFloat, 3},
  F <: AbstractArray{<:AbstractFloat, 5},
  G <: AbstractMatrix{<:AbstractFloat},
  H <: AbstractArray{<:AbstractFloat, 3},
  I <: AbstractArray{<:AbstractFloat, 5},
  J <: AbstractArray{<:AbstractFloat, 3},
  K <: AbstractArray{<:AbstractFloat, 3},
  L <: AbstractArray{<:AbstractFloat, 3},
  M <: AbstractVector{<:AbstractFloat},
}

  # MPI variables.
  comm::A
  master::B
  rank::C
  root::C

  # Start indices of the local grid.
  is::C
  js::C

  # Local grid size.
  nx::C
  ny::C
  nz::C

  # Local grid size with hao/ghost cells
  nxx::C
  nyy::C
  nzz::C

  # Source and destination ranks for halos.
  left::C
  right::C
  back::C
  forw::C

  # Auxiliary arrays for halos.
  send_a2_left::D
  send_a2_right::D
  recv_a2_left::D
  recv_a2_right::D
  send_a3_left::E
  send_a3_right::E
  recv_a3_left::E
  recv_a3_right::E
  send_a5_left::F
  send_a5_right::F
  recv_a5_left::F
  recv_a5_right::F
  send_a2_back::G
  send_a2_forw::G
  recv_a2_back::G
  recv_a2_forw::G
  send_a3_back::H
  send_a3_forw::H
  recv_a3_back::H
  recv_a3_forw::H
  send_a5_back::I
  send_a5_forw::I
  recv_a5_back::I
  recv_a5_forw::I

  # Auxiliary arrays for gather & scatter.
  local_array::J
  master_array::K
  global_array::L

  # Auxiliary arrays for horizontal averaging.
  local_sum::M
  global_sum::M
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
  nbproc = MPI.Comm_size(MPI.COMM_WORLD)

  # Set the start indices.
  istart = -nbx + 1 # left start index of global domain
  jstart = -nby + 1 # back start index of global domain

  # Set number of CPUs in each direction
  if nprocx * nprocy != nbproc
    if master
      println(
        "Virtual topology wrong, nprocx * nprocy should be equal to nbProc!",
      )
      println("[nprocx, nprocy] = [", nprocx, ", ", nprocy, "]")
      println("nbProc = ", nbproc)
    end
    exit()
  end

  # Set dimensions and periodicity.
  dims = [nprocx, nprocy]
  periods = [true, true]

  # Create a Cartesian topology.
  comm = MPI.Cart_create(MPI.COMM_WORLD, dims; periodic = periods)
  rank = MPI.Comm_rank(comm)
  coords = MPI.Cart_coords(comm, rank)
  icoord = coords[1]
  jcoord = coords[2]

  # Set up composition along x coordinate.
  if nprocx == 1
    nx = sizex
    is = istart
  else
    nx1 = div(sizex - 1, nprocx) + 1
    in1 = sizex - nprocx * (nx1 - 1)
    if icoord < in1
      nx = nx1
      is = istart + icoord * nx
    else
      nx = nx1 - 1
      is = istart + in1 * nx1 + (icoord - in1) * nx
    end
  end
  nxx = nx + 2 * nbx + 1
  ie = is + nxx - 1

  # Set up composition along x coordinate.
  if nprocy == 1
    ny = sizey
    js = jstart
  else
    ny1 = div(sizey - 1, nprocy) + 1
    jn1 = sizey - nprocy * (ny1 - 1)
    if jcoord < jn1
      ny = ny1
      js = jstart + jcoord * ny
    else
      ny = ny1 - 1
      js = jstart + jn1 * ny1 + (jcoord - jn1) * ny
    end
  end
  nyy = ny + 2 * nby + 1
  je = js + nyy - 1

  # Set up serial z coordinate.
  nz = sizez
  nzz = nz + 2 * nbz + 1

  # Check consistency between domain size and cpus.
  if sizex < nprocx
    if master
      println(
        "Stopping! Do not use more than ",
        sizex,
        " cpus to compute this application",
      )
    end
    exit()
  end
  if sizey < nprocy
    if master
      println(
        "Stopping! Do not use more than ",
        sizey,
        " cpus to compute this application",
      )
    end
    exit()
  end

  # Number of halo cells should not be larger than number of physical cells.
  if nprocx > 1 && nbx > nx
    if master
      println("Error in init_mpi: nprocx > 1 and nbx > nx!")
      exit()
    end
  end
  if nprocy > 1 && nby > ny
    if master
      println("Error in init_mpi: nprocy > 1 and nby > ny!")
      exit()
    end
  end

  # Find the neighbour processors.
  (left, right) = MPI.Cart_shift(comm, 0, 1)
  (back, forw) = MPI.Cart_shift(comm, 1, 1)

  # Initialize auxiliary arrays for halos.
  (send_a2_left, send_a2_right, recv_a2_left, recv_a2_right) =
    (zeros((nbx, ny + 2 * nby + 1)) for i in 1:4)
  (send_a3_left, send_a3_right, recv_a3_left, recv_a3_right) =
    (zeros((nbx, ny + 2 * nby + 1, nz + 2 * nbz + 1)) for i in 1:4)
  (send_a5_left, send_a5_right, recv_a5_left, recv_a5_right) =
    (zeros((nbx, ny + 2 * nby + 1, nz + 2 * nbz + 1, 3, 2)) for i in 1:4)
  (send_a2_back, send_a2_forw, recv_a2_back, recv_a2_forw) =
    (zeros((nx + 2 * nbx + 1, nby)) for i in 1:4)
  (send_a3_back, send_a3_forw, recv_a3_back, recv_a3_forw) =
    (zeros((nx + 2 * nbx + 1, nby, nz + 2 * nbz + 1)) for i in 1:4)
  (send_a5_back, send_a5_forw, recv_a5_back, recv_a5_forw) =
    (zeros((nx + 2 * nbx + 1, nby, nz + 2 * nbz + 1, 3, 2)) for i in 1:4)

  # Initialize auxiliary arrays for gather & scatter.
  local_array = zeros((nx, ny, nz))
  master_array = zeros((sizex * nprocy, ny, nz))
  global_array = zeros((sizex, sizey, sizez))

  # Initialize auxiliary arrays for horizontal averaging.
  (local_sum, global_sum) = (zeros(nz) for i in 1:2)

  # Return Domain instance.
  return Domain(
    comm,
    master,
    rank,
    root,
    is,
    js,
    nx,
    ny,
    nz,
    nxx,
    nyy,
    nzz,
    left,
    right,
    back,
    forw,
    send_a2_left,
    send_a2_right,
    recv_a2_left,
    recv_a2_right,
    send_a3_left,
    send_a3_right,
    recv_a3_left,
    recv_a3_right,
    send_a5_left,
    send_a5_right,
    recv_a5_left,
    recv_a5_right,
    send_a2_back,
    send_a2_forw,
    recv_a2_back,
    recv_a2_forw,
    send_a3_back,
    send_a3_forw,
    recv_a3_back,
    recv_a3_forw,
    send_a5_back,
    send_a5_forw,
    recv_a5_back,
    recv_a5_forw,
    local_array,
    master_array,
    global_array,
    local_sum,
    global_sum,
  )
end
