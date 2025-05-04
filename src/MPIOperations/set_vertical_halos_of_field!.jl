function set_vertical_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
)
    (; nbz) = namelists.domain
    (;
        comm,
        sizezz,
        nzz,
        ko,
        k0,
        k1,
        down,
        up,
        send_f3_down,
        send_f3_up,
        recv_f3_down,
        recv_f3_up,
    ) = domain

    if ko == 0
        for k in 1:nbz
            @views send_f3_up[:, :, k] .= field[:, :, k1 - k + 1]
        end

        MPI.Sendrecv!(send_f3_up, recv_f3_up, comm; dest = up, source = up)

        for k in 1:nbz
            @views field[:, :, k1 + k] .= recv_f3_up[:, :, k]
        end
    elseif ko + nzz == sizezz
        for k in 1:nbz
            @views send_f3_down[:, :, k] .= field[:, :, k0 + k - 1]
        end

        MPI.Sendrecv!(
            send_f3_down,
            recv_f3_down,
            comm;
            dest = down,
            source = down,
        )

        for k in 1:nbz
            @views field[:, :, k0 - k] .= recv_f3_down[:, :, k]
        end
    else
        for k in 1:nbz
            @views send_f3_up[:, :, k] .= field[:, :, k1 - k + 1]
            @views send_f3_down[:, :, k] .= field[:, :, k0 + k - 1]
        end

        MPI.Sendrecv!(send_f3_up, recv_f3_down, comm; dest = up, source = down)

        MPI.Sendrecv!(send_f3_down, recv_f3_up, comm; dest = down, source = up)

        for k in 1:nbz
            @views field[:, :, k0 - k] .= recv_f3_down[:, :, k]
            @views field[:, :, k1 + k] .= recv_f3_up[:, :, k]
        end
    end

    return
end

function set_vertical_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
)
    (; nbz) = namelists.domain
    (;
        comm,
        sizezz,
        nzz,
        ko,
        k0,
        k1,
        down,
        up,
        send_f5_down,
        send_f5_up,
        recv_f5_down,
        recv_f5_up,
    ) = domain

    if ko == 0
        for k in 1:nbz
            @views send_f5_up[:, :, k, :, :] .= field[:, :, k1 - k + 1, :, :]
        end

        MPI.Sendrecv!(send_f5_up, recv_f5_up, comm; dest = up, source = up)

        for k in 1:nbz
            @views field[:, :, k1 + k, :, :] .= recv_f5_up[:, :, k, :, :]
        end
    elseif ko + nzz == sizezz
        for k in 1:nbz
            @views send_f5_down[:, :, k, :, :] .= field[:, :, k0 + k - 1, :, :]
        end

        MPI.Sendrecv!(
            send_f5_down,
            recv_f5_down,
            comm;
            dest = down,
            source = down,
        )

        for k in 1:nbz
            @views field[:, :, k0 - k, :, :] .= recv_f5_down[:, :, k, :, :]
        end
    else
        for k in 1:nbz
            @views send_f5_up[:, :, k, :, :] .= field[:, :, k1 - k + 1, :, :]
            @views send_f5_down[:, :, k, :, :] .= field[:, :, k0 + k - 1, :, :]
        end

        MPI.Sendrecv!(send_f5_up, recv_f5_down, comm; dest = up, source = down)

        MPI.Sendrecv!(send_f5_down, recv_f5_up, comm; dest = down, source = up)

        for k in 1:nbz
            @views field[:, :, k0 - k, :, :] .= recv_f5_down[:, :, k, :, :]
            @views field[:, :, k1 + k, :, :] .= recv_f5_up[:, :, k, :, :]
        end
    end

    return
end
