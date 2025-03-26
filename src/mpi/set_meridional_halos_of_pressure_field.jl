function set_meridional_halos_of_pressure_field!(
  field::AbstractArray{<:AbstractFloat, 3},
  domain::Domain,
)

  # Get all necessary fields.
  (;
    comm,
    ny,
    back,
    forw,
    send_pressure_back,
    send_pressure_forw,
    recv_pressure_back,
    recv_pressure_forw,
  ) = domain

  # Read slice into contiguous array.
  @views send_pressure_back .= field.parent[:, 2, :]
  @views send_pressure_forw .= field.parent[:, ny + 1, :]

  # back -> forw
  source = back
  dest = forw
  tag = 100

  MPI.Sendrecv!(
    send_pressure_forw,
    recv_pressure_back,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # forw -> back
  source = forw
  dest = back
  tag = 100

  MPI.Sendrecv!(
    send_pressure_back,
    recv_pressure_forw,
    comm;
    dest = dest,
    sendtag = tag,
    source = source,
  )

  # Write auxiliary slice to var field.
  field.parent[:, ny + 2, :] .= recv_pressure_forw
  field.parent[:, 1, :] .= recv_pressure_back

  return
end
