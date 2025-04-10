function set_meridional_halo_ray_volumes!(state::State)
    (; comm, nx, nz, i0, i1, j0, j1, k0, k1, back, forw) = state.domain
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

    @views nray_max_back = maximum(nray[:, j0, :])
    @views nray_max_forw = maximum(nray[:, j1, :])

    nray_max_back = MPI.Allreduce(nray_max_back, +, comm)
    nray_max_forw = MPI.Allreduce(nray_max_forw, +, comm)

    send_rays_forw = zeros(length(fields), nray_max_forw, nx + 2, nz + 2)
    send_rays_back = zeros(length(fields), nray_max_back, nx + 2, nz + 2)

    recv_rays_back = zeros(length(fields), nray_max_forw, nx + 2, nz + 2)
    recv_rays_forw = zeros(length(fields), nray_max_back, nx + 2, nz + 2)

    for (index, field) in enumerate(fields)
        @views send_rays_forw[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j1,
            (k0 - 1):(k1 + 1),
        ]
        @views send_rays_back[index, :, :, :] .= getfield(rays, field)[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j0,
            (k0 - 1):(k1 + 1),
        ]
    end

    MPI.Sendrecv!(
        send_rays_forw,
        recv_rays_back,
        comm;
        dest = forw,
        source = back,
    )

    MPI.Sendrecv!(
        send_rays_back,
        recv_rays_forw,
        comm;
        dest = back,
        source = forw,
    )

    for (index, field) in enumerate(fields)
        @views getfield(rays, field)[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ] = recv_rays_back[index, :, :, :]
        @views getfield(rays, field)[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ] = recv_rays_forw[index, :, :, :]
    end

    @views rays.area_xk[
        1:nray_max_forw,
        (i0 - 1):(i1 + 1),
        j0 - 1,
        (k0 - 1):(k1 + 1),
    ] =
        rays.dxray[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ] * rays.dkray[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_yl[
        1:nray_max_forw,
        (i0 - 1):(i1 + 1),
        j0 - 1,
        (k0 - 1):(k1 + 1),
    ] =
        rays.dyray[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ] * rays.dlray[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_zm[
        1:nray_max_forw,
        (i0 - 1):(i1 + 1),
        j0 - 1,
        (k0 - 1):(k1 + 1),
    ] =
        rays.dzray[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ] * rays.dmray[
            1:nray_max_forw,
            (i0 - 1):(i1 + 1),
            j0 - 1,
            (k0 - 1):(k1 + 1),
        ]

    @views rays.area_xk[
        1:nray_max_back,
        (i0 - 1):(i1 + 1),
        j1 + 1,
        (k0 - 1):(k1 + 1),
    ] =
        rays.dxray[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ] * rays.dkray[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_yl[
        1:nray_max_back,
        (i0 - 1):(i1 + 1),
        j1 + 1,
        (k0 - 1):(k1 + 1),
    ] =
        rays.dyray[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ] * rays.dlray[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ]
    @views rays.area_zm[
        1:nray_max_back,
        (i0 - 1):(i1 + 1),
        j1 + 1,
        (k0 - 1):(k1 + 1),
    ] =
        rays.dzray[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ] * rays.dmray[
            1:nray_max_back,
            (i0 - 1):(i1 + 1),
            j1 + 1,
            (k0 - 1):(k1 + 1),
        ]

    return
end
