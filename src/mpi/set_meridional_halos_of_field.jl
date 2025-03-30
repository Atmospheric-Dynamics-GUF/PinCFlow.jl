function set_meridional_halos_of_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nby) = namelists.domain
  (;
    comm,
    j0,
    j1,
    back,
    forw,
    send_a3_back,
    send_a3_forw,
    recv_a3_back,
    recv_a3_forw,
  ) = domain

  # Read slice into auxiliary array.
  for j in 1:nby
    @views send_a3_forw[:, j, :] .= field[:, j1 - j + 1, :]
    @views send_a3_back[:, j, :] .= field[:, j0 + j - 1, :]
  end

  MPI.Sendrecv!(
    send_a3_forw,
    recv_a3_back,
    comm;
    dest = forw,
    source = back,
  )

  MPI.Sendrecv!(
    send_a3_back,
    recv_a3_forw,
    comm;
    dest = back,
    source = forw,
  )

  # Write auxiliary slice to field.
  for j in 1:nby
    @views field[:, j0 - j, :] .= recv_a3_back[:, j, :]
    @views field[:, j1 + j, :] .= recv_a3_forw[:, j, :]
  end

  return
end

function set_meridional_halos_of_field!(
  field::AbstractArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nby) = namelists.domain
  (;
    comm,
    j0,
    j1,
    back,
    forw,
    send_a5_back,
    send_a5_forw,
    recv_a5_back,
    recv_a5_forw,
  ) = domain

  # Read slice into auxiliary array.
  for j in 1:nby
    @views send_a5_forw[:, j, :, :, :] .= field[:, j1 - j + 1, :, :, :]
    @views send_a5_back[:, j, :, :, :] .= field[:, j0 + j - 1, :, :, :]
  end

  MPI.Sendrecv!(
    send_a5_forw,
    recv_a5_back,
    comm;
    dest = forw,
    source = back,
  )

  MPI.Sendrecv!(
    send_a5_back,
    recv_a5_forw,
    comm;
    dest = back,
    source = forw,
  )

  # Write auxiliary slice to field.
  for j in 1:nby
    @views field[:, j0 - j, :, :, :] .= recv_a5_back[:, j, :, :, :]
    @views field[:, j1 + j, :, :, :] .= recv_a5_forw[:, j, :, :, :]
  end

  return
end
