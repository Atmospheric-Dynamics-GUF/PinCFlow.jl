struct Domain{A <: MPI.Comm, B <: Bool, C <: Integer}

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
  nprocx = nprocx
  nprocy = nprocy
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
    if icoord > in1
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
    if jcoord > jn1
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

  # Return Domain instance.
  return Domain(comm, master, rank, root, is, js, nx, ny, nz, nxx, nyy, nzz)
end
