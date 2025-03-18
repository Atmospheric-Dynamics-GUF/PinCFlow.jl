function set_meridional_halos_of_field!(
  field::OffsetMatrix{<:AbstractFloat},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby) = namelists.domain
  (; comm, nx, ny) = domain

  # Find the neighbour processors and initialize send and reveive buffers.
  (back, forw) = MPI.Cart_shift(comm, 1, 1)
  (ysliceback_send, ysliceforw_send, ysliceback_recv, ysliceforw_recv) = (
    OffsetArray(zeros((nx + 2 * nbx + 1, nby)), (-nbx):(nx + nbx), 1:nby) for
    i in 1:4
  )

  # Set slice size.
  sendcount = nby * (nx + 2 * nbx + 1)
  recvcount = sendcount

  # Read slice into contiguous array.
  for j in 1:nby
    ysliceback_send[:, j] = field[:, j]
    ysliceforw_send[:, j] = field[:, ny - nby + j]
  end

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(ysliceforw_send, ysliceback_recv, comm; dest = dest, sendtag = tag, source = source)

  # forw -> back
  source = forw
  dest = back
  tag = 100

  Mpi.Sendrecv!(ysliceback_send, ysliceforw_recv, comm; dest = dest, sendtag = tag, source = source)

  # write auxiliary slice to var field
  for j in 1:nby
    field[:, ny + j] = ysliceforw_recv[:, j]
    field[:, -nby + j] = ysliceback_recv[:, j]
  end

  return
end

function set_meridional_halos_of_field!(
  field::OffsetArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (; comm, nx, ny, nz) = domain

  # Find the neighbour processors and initialize send and reveive buffers.
  (back, forw) = MPI.Cart_shift(comm, 1, 1)
  (ysliceback_send, ysliceforw_send, ysliceback_recv, ysliceforw_recv) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, nby, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      1:nby,
      (-nbz):(nz + nbz),
    ) for i in 1:4
  )

  # Set slice size.
  sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1)
  recvcount = sendcount

  # Read slice into contiguous array.
  for j in 1:nby
    ysliceback_send[:, j, :] = field[:, j, :]
    ysliceforw_send[:, j, :] = field[:, ny - nby + j, :]
  end

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(ysliceforw_send, ysliceback_recv, comm; dest = dest, sendtag = tag, source = source)

  # forw -> back
  source = forw
  dest = back
  tag = 100

  Mpi.Sendrecv!(ysliceback_send, ysliceforw_recv, comm; dest = dest, sendtag = tag, source = source)

  # write auxiliary slice to var field
  for j in 1:nby
    field[:, ny + j, :] = ysliceforw_recv[:, j, :]
    field[:, -nby + j, :] = ysliceback_recv[:, j, :]
  end

  return
end

function set_meridional_halos_of_field!(
  field::OffsetArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)
  (back, forw) = MPI.Cart_shift(comm, 1, 1)
  (ysliceback_send, ysliceforw_send, ysliceback_recv, ysliceforw_recv) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, nby, nz + 2 * nbz + 1, 3, 2)),
      (-nbx):(nx + nbx),
      1:nby,
      (-nbz):(nz + nbz),
      1:3,
      0:1,
    ) for i in 1:4
  )

  # Set slice size.
  sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1) * 6
  recvcount = sendcount

  # Read slice into contiguous array.
  for j in 1:nby
    ysliceback_send[:, j, :, :, :] = field[:, j, :, :, :]
    ysliceforw_send[:, j, :, :, :] = field[:, ny - nby + j, :, :, :]
  end

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(ysliceforw_send, ysliceback_recv, comm; dest = dest, sendtag = tag, source = source)

  # forw -> back
  source = forw
  dest = back
  tag = 100

  Mpi.Sendrecv!(ysliceback_send, ysliceforw_recv, comm; dest = dest, sendtag = tag, source = source)

  # write auxiliary slice to var field
  for j in 1:nby
    field[:, ny + j, :, :, :] = ysliceforw_recv[:, j, :, :, :]
    field[:, -nby + j, :, :, :] = ysliceback_recv[:, j, :, :, :]
  end

  return
end
