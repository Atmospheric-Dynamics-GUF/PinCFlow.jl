function set_zonal_halos_of_field!(
  field::OffsetMatrix{<:AbstractFloat},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby) = namelists.domain
  (; comm, nx, ny) = domain

  # Find the neighbour processors and initialize send and reveive buffers.
  (left, right) = MPI.Cart_shift(comm, 0, 1)
  (xsliceleft_send, xsliceright_send, xsliceleft_recv, xsliceright_recv) =
    (zeros((nbx, ny + 2 * nby + 1)) for i in 1:4)

  # Set slice size.
  sendcount = nbx * (ny + 2 * nby + 1)
  recvcount = sendcount

  # Read slice into contiguous array
  for i in 1:nbx
    xsliceleft_send[i, :] = field[i, :]
    xsliceright_send[i, :] = field[nx - nbx + i, :]
  end

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    xsliceright_send,
    xsliceleft_recv,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # right -> left
  source = right
  dest = left
  tag = 100

  MPI.Sendrecv!(
    xsliceleft_send,
    xsliceright_recv,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for i in 1:nbx
    field[nx + i, :] = xsliceright_recv[i, :]
    field[-nbx + i, :] = xsliceleft_recv[i, :]
  end

  return
end

function set_zonal_halos_of_field!(
  field::OffsetArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (; comm, nx, ny, nz) = domain

  # Find the neighbour processors and initialize send and reveive buffers.
  (left, right) = MPI.Cart_shift(comm, 0, 1)
  (xsliceleft_send, xsliceright_send, xsliceleft_recv, xsliceright_recv) =
    (zeros((nbx, ny + 2 * nby + 1, nz + 2 * nbz + 1)) for i in 1:4)

  # Set slice size.
  sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1)
  recvcount = sendcount

  # Read slice into contiguous array
  for i in 1:nbx
    xsliceleft_send[i, :, :] = field[i, :, :]
    xsliceright_send[i, :, :] = field[nx - nbx + i, :, :]
  end

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    xsliceright_send,
    xsliceleft_recv,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # right -> left
  source = right
  dest = left
  tag = 100

  MPI.Sendrecv!(
    xsliceleft_send,
    xsliceright_recv,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for i in 1:nbx
    field[nx + i, :, :] = xsliceright_recv[i, :, :]
    field[-nbx + i, :, :] = xsliceleft_recv[i, :, :]
  end

  return
end

function set_zonal_halos_of_field!(
  field::OffsetArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (; comm, nx, ny, nz) = domain

  # Find the neighbour processors and initialize send and reveive buffers.
  (left, right) = MPI.Cart_shift(comm, 0, 1)
  (xsliceleft_send, xsliceright_send, xsliceleft_recv, xsliceright_recv) =
    (zeros((nbx, ny + 2 * nby + 1, nz + 2 * nbz + 1, 3, 2)) for i in 1:4)

  # Set slice size.
  sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1) * 6
  recvcount = sendcount

  # Read slice into contiguous array
  for i in 1:nbx
    xsliceleft_send[i, :, :, :, :] = field[i, :, :, :, :]
    xsliceright_send[i, :, :, :, :] = field[nx - nbx + i, :, :, :, :]
  end

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    xsliceright_send,
    xsliceleft_recv,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # right -> left
  source = right
  dest = left
  tag = 100

  MPI.Sendrecv!(
    xsliceleft_send,
    xsliceright_recv,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for i in 1:nbx
    field[nx + i, :, :, :, :] = xsliceright_recv[i, :, :, :, :]
    field[-nbx + i, :, :, :, :] = xsliceleft_recv[i, :, :, :, :]
  end

  return
end
