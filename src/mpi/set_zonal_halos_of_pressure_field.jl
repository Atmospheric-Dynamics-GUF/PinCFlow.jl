function set_zonal_halos_of_pressure_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  domain::Domain,
)

  # Get all necessary fields.
  (;
    comm,
    nx,
    left,
    right,
    send_pressure_left,
    send_pressure_right,
    recv_pressure_left,
    recv_pressure_right,
  ) = domain

  # Read slice into contiguous array
  @views send_pressure_left .= field.parent[2, :, :]
  @views send_pressure_right .= field.parent[nx + 1, :, :]

  # left -> right
  source = left
  dest = right
  tag = 100

  MPI.Sendrecv!(
    send_pressure_right,
    recv_pressure_left,
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
    send_pressure_left,
    recv_pressure_right,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # Write auxiliary slice to var field.
  field.parent[nx + 2, :, :] .= recv_pressure_right
  field.parent[1, :, :] .= recv_pressure_left

  return
end
