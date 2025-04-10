function set_zonal_halo_ray_volumes!(state::State)
    (; comm, ny, nz, i0, i1, j0, j1, k0, k1, left, right) = state.domain
    (; nray, rays) = state.wkb

    fields = (
        :x,
        :y,
        :z,
        :k,
        :l,
        :m,
        :dxray,
        :dyray,
        :dzray,
        :dkray,
        :dlray,
        :dmray,
        :omega,
        :dens,
    )

    @views nray_max_left = maximum(nray[i0, :, :])
    @views nray_max_right = maximum(nray[i1, :, :])

    nray_max_left = MPI.Allreduce(nray_max_left, +, comm)
    nray_max_right = MPI.Allreduce(nray_max_right, +, comm)

    send_rays_right = zeros(length(fields), nray_max_right, ny + 2, nz + 2)
    send_rays_left = zeros(length(fields), nray_max_left, ny + 2, nz + 2)

    recv_rays_left = zeros(length(fields), nray_max_right, ny + 2, nz + 2)
    recv_rays_right = zeros(length(fields), nray_max_left, ny + 2, nz + 2)

    for (index, field) in enumerate(fields)
        @views send_rays_right[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_right,
            i1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
        @views send_rays_left[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_left,
            i0,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
    end

    MPI.Sendrecv!(
        send_rays_right,
        recv_rays_left,
        comm;
        dest = right,
        source = left,
    )

    MPI.Sendrecv!(
        send_rays_left,
        recv_rays_right,
        comm;
        dest = left,
        source = right,
    )

    for (index, field) in enumerate(fields)
        @views getfield(rays, field)[
            1:nray_max_left,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] = recv_rays_left[index, :, :, :]
        @views getfield(rays, field)[
            1:nray_max_right,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] = recv_rays_right[index, :, :, :]
    end

    @views rays.area_xk[
        1:nray_max_right,
        i0 - 1,
        (j0 - 1):(j1 + 1),
        (k0 - 1):(k1 + 1),
    ] =
        rays.dxray[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] * rays.dkray[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_yl[
        1:nray_max_right,
        i0 - 1,
        (j0 - 1):(j1 + 1),
        (k0 - 1):(k1 + 1),
    ] =
        rays.dyray[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] * rays.dlray[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_zm[
        1:nray_max_right,
        i0 - 1,
        (j0 - 1):(j1 + 1),
        (k0 - 1):(k1 + 1),
    ] =
        rays.dzray[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] * rays.dmray[
            1:nray_max_right,
            i0 - 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]

    @views rays.area_xk[
        1:nray_max_left,
        i1 + 1,
        (j0 - 1):(j1 + 1),
        (k0 - 1):(k1 + 1),
    ] =
        rays.dxray[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] * rays.dkray[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_yl[
        1:nray_max_left,
        i1 + 1,
        (j0 - 1):(j1 + 1),
        (k0 - 1):(k1 + 1),
    ] =
        rays.dyray[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] * rays.dlray[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_zm[
        1:nray_max_left,
        i1 + 1,
        (j0 - 1):(j1 + 1),
        (k0 - 1):(k1 + 1),
    ] =
        rays.dzray[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ] * rays.dmray[
            1:nray_max_left,
            i1 + 1,
            (j0 - 1):(j1 + 1),
            (k0 - 1):(k1 + 1),
        ]

    return
end
