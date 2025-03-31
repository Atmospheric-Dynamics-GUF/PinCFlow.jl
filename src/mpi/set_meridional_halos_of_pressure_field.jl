function set_meridional_halos_of_pressure_field!(
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
    back,
    forw,
    send_pressure_back,
    send_pressure_forw,
    recv_pressure_back,
    recv_pressure_forw,
  ) = domain

  # Read slice into auxiliary array.
  @views send_pressure_forw .= field[(i0 - 1):(i1 + 1), j1, (k0 - 1):(k1 + 1)]
  @views send_pressure_back .= field[(i0 - 1):(i1 + 1), j0, (k0 - 1):(k1 + 1)]

  MPI.Sendrecv!(
    send_pressure_forw,
    recv_pressure_back,
    comm;
    dest = forw,
    source = back,
  )

  MPI.Sendrecv!(
    send_pressure_back,
    recv_pressure_forw,
    comm;
    dest = back,
    source = forw,
  )

  # Write auxiliary slice to field.
  field[(i0 - 1):(i1 + 1), j0 - 1, (k0 - 1):(k1 + 1)] .= recv_pressure_back
  field[(i0 - 1):(i1 + 1), j1 + 1, (k0 - 1):(k1 + 1)] .= recv_pressure_forw

  return
end
