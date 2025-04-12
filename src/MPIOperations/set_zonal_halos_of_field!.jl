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
        send_f3_left,
        send_f3_right,
        recv_f3_left,
        recv_f3_right,
    ) = domain

    # Read slice into auxiliary array.
    for i in 1:nbx
        @views send_f3_right[i, :, :] .= field[i1 - i + 1, :, :]
        @views send_f3_left[i, :, :] .= field[i0 + i - 1, :, :]
    end

    MPI.Sendrecv!(
        send_f3_right,
        recv_f3_left,
        comm;
        dest = right,
        source = left,
    )

    MPI.Sendrecv!(
        send_f3_left,
        recv_f3_right,
        comm;
        dest = left,
        source = right,
    )

    # Write auxiliary slice to field.
    for i in 1:nbx
        @views field[i0 - i, :, :] .= recv_f3_left[i, :, :]
        @views field[i1 + i, :, :] .= recv_f3_right[i, :, :]
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
        send_f5_left,
        send_f5_right,
        recv_f5_left,
        recv_f5_right,
    ) = domain

    # Read slice into auxiliary array.
    for i in 1:nbx
        @views send_f5_right[i, :, :, :, :] .= field[i1 - i + 1, :, :, :, :]
        @views send_f5_left[i, :, :, :, :] .= field[i0 + i - 1, :, :, :, :]
    end

    MPI.Sendrecv!(
        send_f5_right,
        recv_f5_left,
        comm;
        dest = right,
        source = left,
    )

    MPI.Sendrecv!(
        send_f5_left,
        recv_f5_right,
        comm;
        dest = left,
        source = right,
    )

    # Write auxiliary slice to field.
    for i in 1:nbx
        @views field[i0 - i, :, :, :, :] .= recv_f5_left[i, :, :, :, :]
        @views field[i1 + i, :, :, :, :] .= recv_f5_right[i, :, :, :, :]
    end

    return
end
