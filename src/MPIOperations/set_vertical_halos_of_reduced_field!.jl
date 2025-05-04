function set_vertical_halos_of_reduced_field!(
    field::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
)

    # Get all necessary fields.
    (;
        comm,
        sizezz,
        nzz,
        ko,
        k0,
        k1,
        down,
        up,
        send_rf3_down,
        send_rf3_up,
        recv_rf3_down,
        recv_rf3_up,
    ) = domain

    if ko == 0
        @views send_rf3_up .= field[:, :, k1]

        # MPI.Send(send_rf3_up, comm; dest = up)
        # MPI.Recv!(recv_rf3_up, comm; source = up)
        MPI.Sendrecv!(send_rf3_up, recv_rf3_up, comm; dest = up, source = up)

        field[:, :, k1 + 1] .= recv_rf3_up
    elseif ko + nzz == sizezz
        @views send_rf3_down .= field[:, :, k0]

        # MPI.Send(send_rf3_down, comm; dest = down)
        # MPI.Recv!(recv_rf3_down, comm; source = down)
        MPI.Sendrecv!(
            send_rf3_down,
            recv_rf3_down,
            comm;
            dest = down,
            source = down,
        )

        field[:, :, k0 - 1] .= recv_rf3_down
    else
        @views send_rf3_up .= field[:, :, k1]
        @views send_rf3_down .= field[:, :, k0]

        MPI.Sendrecv!(
            send_rf3_up,
            recv_rf3_down,
            comm;
            dest = up,
            source = down,
        )

        MPI.Sendrecv!(
            send_rf3_down,
            recv_rf3_up,
            comm;
            dest = down,
            source = up,
        )

        field[:, :, k0 - 1] .= recv_rf3_down
        field[:, :, k1 + 1] .= recv_rf3_up
    end

    return
end
