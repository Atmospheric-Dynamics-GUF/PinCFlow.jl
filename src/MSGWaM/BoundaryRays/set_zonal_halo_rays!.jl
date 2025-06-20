"""
    set_zonal_halo_rays!(state::State)

    Synchronize ray volumes across zonal domain boundaries. 

# Arguments

  - `state::State`: Complete simulation state
"""
function set_zonal_halo_rays!(state::State)
    (; comm, ny, nz, i0, i1, j0, j1, k0, k1, left, right) = state.domain
    (; nray, rays) = state.wkb

    fields = fieldnames(Rays)

    @views nray_max_left =
        maximum(nray[i0, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)])
    @views nray_max_right =
        maximum(nray[i1, (j0 - 1):(j1 + 1), (k0 - 1):(k1 + 1)])

    nray_max_left = MPI.Allreduce(nray_max_left, max, comm)
    nray_max_right = MPI.Allreduce(nray_max_right, max, comm)

    send_right = zeros(length(fields), nray_max_right, ny + 2, nz + 2)
    send_left = zeros(length(fields), nray_max_left, ny + 2, nz + 2)

    receive_left = zeros(length(fields), nray_max_right, ny + 2, nz + 2)
    receive_right = zeros(length(fields), nray_max_left, ny + 2, nz + 2)

    for (index, field) in enumerate(fields)
        @views send_right[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_right,
            i1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
        @views send_left[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_left,
            i0,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
    end

    MPI.Sendrecv!(send_right, receive_left, comm; dest = right, source = left)

    MPI.Sendrecv!(send_left, receive_right, comm; dest = left, source = right)

    for (index, field) in enumerate(fields)
        @views getfield(rays, field)[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] .= receive_left[index, :, :, :]
        @views getfield(rays, field)[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] .= receive_right[index, :, :, :]
    end

    return
end
