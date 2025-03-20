function set_zonal_halos_of_field!(
  field::AbstractMatrix{<:AbstractFloat},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby) = namelists.domain
  (;
    comm,
    nx,
    ny,
    left,
    right,
    send_a2_left,
    send_a2_right,
    recv_a2_left,
    recv_a2_right,
  ) = domain

  # Read slice into contiguous array
  for i in 1:nbx
    send_a2_left[i, :] = field[i, :]
    send_a2_right[i, :] = field[nx - nbx + i, :]
  end

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    send_a2_right,
    recv_a2_left,
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
    send_a2_left,
    recv_a2_right,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for i in 1:nbx
    field[nx + i, :] = recv_a2_right[i, :]
    field[-nbx + i, :] = recv_a2_left[i, :]
  end

  return
end

function set_zonal_halos_of_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (;
    comm,
    nx,
    ny,
    nz,
    left,
    right,
    send_a3_left,
    send_a3_right,
    recv_a3_left,
    recv_a3_right,
  ) = domain

  # Read slice into contiguous array
  for i in 1:nbx
    send_a3_left[i, :, :] = field[i, :, :]
    send_a3_right[i, :, :] = field[nx - nbx + i, :, :]
  end

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    send_a3_right,
    recv_a3_left,
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
    send_a3_left,
    recv_a3_right,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for i in 1:nbx
    field[nx + i, :, :] = recv_a3_right[i, :, :]
    field[-nbx + i, :, :] = recv_a3_left[i, :, :]
  end

  return
end

function set_zonal_halos_of_field!(
  field::AbstractArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx, nby, nbz) = namelists.domain
  (;
    comm,
    nx,
    ny,
    nz,
    left,
    right,
    send_a5_left,
    send_a5_right,
    recv_a5_left,
    recv_a5_right,
  ) = domain

  # Read slice into contiguous array
  for i in 1:nbx
    send_a5_left[i, :, :, :, :] = field[i, :, :, :, :]
    send_a5_right[i, :, :, :, :] = field[nx - nbx + i, :, :, :, :]
  end

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    send_a5_right,
    recv_a5_left,
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
    send_a5_left,
    recv_a5_right,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for i in 1:nbx
    field[nx + i, :, :, :, :] = recv_a5_right[i, :, :, :, :]
    field[-nbx + i, :, :, :, :] = recv_a5_left[i, :, :, :, :]
  end

  return
end
