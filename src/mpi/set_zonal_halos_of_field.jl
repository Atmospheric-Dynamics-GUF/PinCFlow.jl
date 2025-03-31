function set_zonal_halos_of_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx) = namelists.domain
  (;
    comm,
    i0,
    i1,
    left,
    right,
    send_a3_left,
    send_a3_right,
    recv_a3_left,
    recv_a3_right,
  ) = domain

  # Read slice into auxiliary array.
  for i in 1:nbx
    @views send_a3_right[i, :, :] .= field[i1 - i + 1, :, :]
    @views send_a3_left[i, :, :] .= field[i0 + i - 1, :, :]
  end

  MPI.Sendrecv!(send_a3_right, recv_a3_left, comm; dest = right, source = left)

  MPI.Sendrecv!(send_a3_left, recv_a3_right, comm; dest = left, source = right)

  # Write auxiliary slice to field.
  for i in 1:nbx
    @views field[i0 - i, :, :] .= recv_a3_left[i, :, :]
    @views field[i1 + i, :, :] .= recv_a3_right[i, :, :]
  end

  return
end

function set_zonal_halos_of_field!(
  field::AbstractArray{<:AbstractFloat, 5},
  namelists::Namelists,
  domain::Domain,
)

  # Get all necessary fields.
  (; nbx) = namelists.domain
  (;
    comm,
    i0,
    i1,
    left,
    right,
    send_a5_left,
    send_a5_right,
    recv_a5_left,
    recv_a5_right,
  ) = domain

  # Read slice into auxiliary array.
  for i in 1:nbx
    @views send_a5_right[i, :, :, :, :] .= field[i1 - i + 1, :, :, :, :]
    @views send_a5_left[i, :, :, :, :] .= field[i0 + i - 1, :, :, :, :]
  end

  MPI.Sendrecv!(send_a5_right, recv_a5_left, comm; dest = right, source = left)

  MPI.Sendrecv!(send_a5_left, recv_a5_right, comm; dest = left, source = right)

  # Write auxiliary slice to field.
  for i in 1:nbx
    @views field[i0 - i, :, :, :, :] .= recv_a5_left[i, :, :, :, :]
    @views field[i1 + i, :, :, :, :] .= recv_a5_right[i, :, :, :, :]
  end

  return
end
