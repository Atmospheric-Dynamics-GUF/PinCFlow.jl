function set_vertical_halos_of_reduced_field!(
    field::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
)
    (;
        comm,
        sizezz,
        nzz,
        ko,
        i0,
        i1,
        j0,
        j1,
        k0,
        k1,
        down,
        up,
        send_rf3_down,
        send_rf3_up,
        recv_rf3_down,
        recv_rf3_up,
    ) = domain

    ix = (i0 - 1):(i1 + 1)
    jy = (j0 - 1):(j1 + 1)

    if ko == 0
        @views send_rf3_up .= field[ix, jy, k1]

        MPI.Sendrecv!(send_rf3_up, recv_rf3_up, comm; dest = up, source = up)

        field[ix, jy, k1 + 1] .= recv_rf3_up
    elseif ko + nzz == sizezz
        @views send_rf3_down .= field[ix, jy, k0]

        MPI.Sendrecv!(
            send_rf3_down,
            recv_rf3_down,
            comm;
            dest = down,
            source = down,
        )

        field[ix, jy, k0 - 1] .= recv_rf3_down
    else
        @views send_rf3_up .= field[ix, jy, k1]
        @views send_rf3_down .= field[ix, jy, k0]

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

        field[ix, jy, k0 - 1] .= recv_rf3_down
        field[ix, jy, k1 + 1] .= recv_rf3_up
    end

    return
end
