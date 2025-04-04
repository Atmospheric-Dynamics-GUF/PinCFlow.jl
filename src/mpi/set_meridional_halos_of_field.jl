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
    send_f3_back,
    send_f3_forw,
    recv_f3_back,
    recv_f3_forw,
  ) = domain

  # Read slice into auxiliary array.
  for j in 1:nby
    @views send_f3_forw[:, j, :] .= field[:, j1 - j + 1, :]
    @views send_f3_back[:, j, :] .= field[:, j0 + j - 1, :]
  end

  MPI.Sendrecv!(send_f3_forw, recv_f3_back, comm; dest = forw, source = back)

  MPI.Sendrecv!(send_f3_back, recv_f3_forw, comm; dest = back, source = forw)

  # Write auxiliary slice to field.
  for j in 1:nby
    @views field[:, j0 - j, :] .= recv_f3_back[:, j, :]
    @views field[:, j1 + j, :] .= recv_f3_forw[:, j, :]
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
    send_f5_back,
    send_f5_forw,
    recv_f5_back,
    recv_f5_forw,
  ) = domain

  # Read slice into auxiliary array.
  for j in 1:nby
    @views send_f5_forw[:, j, :, :, :] .= field[:, j1 - j + 1, :, :, :]
    @views send_f5_back[:, j, :, :, :] .= field[:, j0 + j - 1, :, :, :]
  end

  MPI.Sendrecv!(send_f5_forw, recv_f5_back, comm; dest = forw, source = back)

  MPI.Sendrecv!(send_f5_back, recv_f5_forw, comm; dest = back, source = forw)

  # Write auxiliary slice to field.
  for j in 1:nby
    @views field[:, j0 - j, :, :, :] .= recv_f5_back[:, j, :, :, :]
    @views field[:, j1 + j, :, :, :] .= recv_f5_forw[:, j, :, :, :]
  end

  return
end
