function set_meridional_halos_of_field!(
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
    back,
    forw,
    send_a2_back,
    send_a2_forw,
    recv_a2_back,
    recv_a2_forw,
  ) = domain

  # Read slice into contiguous array.
  for j in 1:nby
    @views send_a2_back[:, j] = field[:, j]
    @views send_a2_forw[:, j] = field[:, ny - nby + j]
  end

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(
    send_a2_forw,
    recv_a2_back,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # forw -> back
  source = forw
  dest = back
  tag = 100

  Mpi.Sendrecv!(
    send_a2_back,
    recv_a2_forw,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for j in 1:nby
    @views field[:, ny + j] = recv_a2_forw[:, j]
    @views field[:, -nby + j] = recv_a2_back[:, j]
  end

  return
end

function set_meridional_halos_of_field!(
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
    back,
    forw,
    send_a3_back,
    send_a3_forw,
    recv_a3_back,
    recv_a3_forw,
  ) = domain

  # Read slice into contiguous array.
  for j in 1:nby
    @views send_a3_back[:, j, :] = field[:, j, :]
    @views send_a3_forw[:, j, :] = field[:, ny - nby + j, :]
  end

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(
    send_a3_forw,
    recv_a3_back,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # forw -> back
  source = forw
  dest = back
  tag = 100

  Mpi.Sendrecv!(
    send_a3_back,
    recv_a3_forw,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for j in 1:nby
    @views field[:, ny + j, :] = recv_a3_forw[:, j, :]
    @views field[:, -nby + j, :] = recv_a3_back[:, j, :]
  end

  return
end

function set_meridional_halos_of_field!(
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
    back,
    forw,
    send_a5_back,
    send_a5_forw,
    recv_a5_back,
    recv_a5_forw,
  ) = domain

  # Read slice into contiguous array.
  for j in 1:nby
    @views send_a5_back[:, j, :, :, :] = field[:, j, :, :, :]
    @views send_a5_forw[:, j, :, :, :] = field[:, ny - nby + j, :, :, :]
  end

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(
    send_a5_forw,
    recv_a5_back,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # forw -> back
  source = forw
  dest = back
  tag = 100

  Mpi.Sendrecv!(
    send_a5_back,
    recv_a5_forw,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # write auxiliary slice to var field
  for j in 1:nby
    @views field[:, ny + j, :, :, :] = recv_a5_forw[:, j, :, :, :]
    @views field[:, -nby + j, :, :, :] = recv_a5_back[:, j, :, :, :]
  end

  return
end
