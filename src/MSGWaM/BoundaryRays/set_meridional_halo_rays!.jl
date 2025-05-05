function set_meridional_halo_rays!(state::State)
    (; comm, nx, nz, i0, i1, j0, j1, k0, k1, back, forw) = state.domain
    (; nray, rays) = state.wkb

    fields = fieldnames(Rays)

    @views nray_max_back =
        maximum(nray[(i0 - 1):(i1 + 1), j0, (k0 - 1):(k1 + 1)])
    @views nray_max_forw =
        maximum(nray[(i0 - 1):(i1 + 1), j1, (k0 - 1):(k1 + 1)])

    nray_max_back = MPI.Allreduce(nray_max_back, max, comm)
    nray_max_forw = MPI.Allreduce(nray_max_forw, max, comm)

    send_forw = zeros(length(fields), nray_max_forw, nx + 2, nz + 2)
    send_back = zeros(length(fields), nray_max_back, nx + 2, nz + 2)

    receive_back = zeros(length(fields), nray_max_forw, nx + 2, nz + 2)
    receive_forw = zeros(length(fields), nray_max_back, nx + 2, nz + 2)

    for (index, field) in enumerate(fields)
        @views send_forw[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j1,
            (k0 - 1):(k1 + 1),
        ]
        @views send_back[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j0,
            (k0 - 1):(k1 + 1),
        ]
    end

    MPI.Sendrecv!(
        send_forw,
        receive_back,
        comm;
        dest = forw,
        source = back,
    )

    MPI.Sendrecv!(
        send_back,
        receive_forw,
        comm;
        dest = back,
        source = forw,
    )

    for (index, field) in enumerate(fields)
        @views getfield(rays, field)[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ] .= receive_back[index, :, :, :]
        @views getfield(rays, field)[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ] .= receive_forw[index, :, :, :]
    end

    return
end
