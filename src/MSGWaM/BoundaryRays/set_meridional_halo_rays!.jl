"""
    set_meridional_halo_rays!(state::State)

# Arguments

  - `state::State`: Complete simulation state
"""
function set_meridional_halo_rays!(state::State)
    (; comm, nx, nz, i0, i1, j0, j1, k0, k1, backward, forward) = state.domain
    (; nray, rays) = state.wkb

    fields = fieldnames(Rays)

    @views nray_max_backward =
        maximum(nray[(i0 - 1):(i1 + 1), j0, (k0 - 1):(k1 + 1)])
    @views nray_max_forward =
        maximum(nray[(i0 - 1):(i1 + 1), j1, (k0 - 1):(k1 + 1)])

    nray_max_backward = MPI.Allreduce(nray_max_backward, max, comm)
    nray_max_forward = MPI.Allreduce(nray_max_forward, max, comm)

    send_forward = zeros(length(fields), nray_max_forward, nx + 2, nz + 2)
    send_backward = zeros(length(fields), nray_max_backward, nx + 2, nz + 2)

    receive_backward = zeros(length(fields), nray_max_forward, nx + 2, nz + 2)
    receive_forward = zeros(length(fields), nray_max_backward, nx + 2, nz + 2)

    for (index, field) in enumerate(fields)
        @views send_forward[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_forward,
            (i0 - 1):(i1 + 1),
            j1,
            (k0 - 1):(k1 + 1),
        ]
        @views send_backward[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_backward,
            (i0 - 1):(i1 + 1),
            j0,
            (k0 - 1):(k1 + 1),
        ]
    end

    MPI.Sendrecv!(
        send_forward,
        receive_backward,
        comm;
        dest = forward,
        source = backward,
    )

    MPI.Sendrecv!(
        send_backward,
        receive_forward,
        comm;
        dest = backward,
        source = forward,
    )

    for (index, field) in enumerate(fields)
        @views getfield(rays, field)[
            1:nray_max_forward,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ] .= receive_backward[index, :, :, :]
        @views getfield(rays, field)[
            1:nray_max_backward,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ] .= receive_forward[index, :, :, :]
    end

    return
end
