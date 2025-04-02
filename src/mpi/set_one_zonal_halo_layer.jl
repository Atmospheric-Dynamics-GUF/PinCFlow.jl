function set_one_zonal_halo_layer!(
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
    send_o3_left,
    send_o3_right,
    recv_o3_left,
    recv_o3_right,
  ) = domain

  # Read slice into auxiliary array
  @views send_o3_right .= field[i1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]
  @views send_o3_left .= field[i0, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)]

  MPI.Sendrecv!(
    send_o3_right,
    recv_o3_left,
    comm;
    dest = right,
    source = left,
  )

  MPI.Sendrecv!(
    send_o3_left,
    recv_o3_right,
    comm;
    dest = left,
    source = right,
  )

  # Write auxiliary slice to var field.
  field[i0 - 1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .= recv_o3_left
  field[i1 + 1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)] .= recv_o3_right

  return
end
