function set_zonal_halos_of_pressure_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  domain::Domain,
)

  # Get all necessary fields.
  (;
    comm,
    i0,
    i1,
    j0,
    j1,
    k0,
    k1,
    left,
    right,
    send_pressure_left,
    send_pressure_right,
    recv_pressure_left,
    recv_pressure_right,
  ) = domain

  # Read slice into auxiliary array
  @views send_pressure_right .= field[i1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]
  @views send_pressure_left .= field[i0, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]

  MPI.Sendrecv!(
    send_pressure_right,
    recv_pressure_left,
    comm;
    dest = right,
    source = left,
  )

  MPI.Sendrecv!(
    send_pressure_left,
    recv_pressure_right,
    comm;
    dest = left,
    source = right,
  )

  # Write auxiliary slice to var field.
  field[i0 - 1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .= recv_pressure_left
  field[i1 + 1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .= recv_pressure_right

  return
end
