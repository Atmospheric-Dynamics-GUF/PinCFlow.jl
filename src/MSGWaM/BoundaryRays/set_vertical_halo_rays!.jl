function set_vertical_halo_rays!(state::State)
    (; comm, sizezz, nzz, nx, ny, ko, i0, i1, j0, j1, k0, k1, down, up) =
        state.domain
    (; nray, rays) = state.wkb

    fields = fieldnames(Rays)

    @views nray_max_down =
        maximum(nray[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), k0])
    @views nray_max_up = maximum(nray[(i0 - 1):(i1 + 1), (j0 - 1):(j1 + 1), k1])

    nray_max_down = MPI.Allreduce(nray_max_down, max, comm)
    nray_max_up = MPI.Allreduce(nray_max_up, max, comm)

    send_up = zeros(length(fields), nray_max_up, nx + 2, ny + 2)
    send_down = zeros(length(fields), nray_max_down, nx + 2, ny + 2)

    receive_down = zeros(length(fields), nray_max_up, nx + 2, ny + 2)
    receive_up = zeros(length(fields), nray_max_down, nx + 2, ny + 2)

    if ko == 0
        for (index, field) in enumerate(fields)
            @views send_up[index, :, :, :] .= getfield(rays, field)[
                1:nray_max_up,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k1,
            ]
        end

        MPI.Sendrecv!(send_up, receive_up, comm; dest = up, source = up)

        for (index, field) in enumerate(fields)
            @views getfield(rays, field)[
                1:nray_max_down,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k1 + 1,
            ] .= receive_up[index, :, :, :]
        end
    elseif ko + nzz == sizezz
        for (index, field) in enumerate(fields)
            @views send_down[index, :, :, :] .= getfield(rays, field)[
                1:nray_max_down,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k0,
            ]
        end

        MPI.Sendrecv!(
            send_down,
            receive_down,
            comm;
            dest = down,
            source = down,
        )

        for (index, field) in enumerate(fields)
            @views getfield(rays, field)[
                1:nray_max_up,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k0 - 1,
            ] .= receive_down[index, :, :, :]
        end
    else
        for (index, field) in enumerate(fields)
            @views send_up[index, :, :, :] .= getfield(rays, field)[
                1:nray_max_up,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k1,
            ]
            @views send_down[index, :, :, :] .= getfield(rays, field)[
                1:nray_max_down,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k0,
            ]
        end

        MPI.Sendrecv!(
            send_up,
            receive_down,
            comm;
            dest = up,
            source = down,
        )

        MPI.Sendrecv!(
            send_down,
            receive_up,
            comm;
            dest = down,
            source = up,
        )

        for (index, field) in enumerate(fields)
            @views getfield(rays, field)[
                1:nray_max_up,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k0 - 1,
            ] .= receive_down[index, :, :, :]
            @views getfield(rays, field)[
                1:nray_max_down,
                (i0 - 1):(i1 + 1),
                (j0 - 1):(j1 + 1),
                k1 + 1,
            ] .= receive_up[index, :, :, :]
        end
    end

    return
end
