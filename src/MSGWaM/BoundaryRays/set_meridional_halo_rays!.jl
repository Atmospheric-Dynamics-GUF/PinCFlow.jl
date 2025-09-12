"""
```julia
set_meridional_halo_rays!(state::State)
```

Exchange ray volumes in meridional halo cells.

Performs bidirectional MPI communication between backward and forward neighbor processes. The number of communicated ray volumes is determined from the maximum counts at the backward and forward boundaries of the MPI subdomains.

# Arguments

  - `state`: Model state.
"""
function set_meridional_halo_rays! end

function set_meridional_halo_rays!(state::State)
    (; comm, nx, nz, i0, i1, j0, j1, k0, k1, backward, forward) = state.domain
    (; nray, rays) = state.wkb

    fields = fieldnames(Rays)

    ii = (i0 - 1):(i1 + 1)
    kk = (k0 - 1):(k1 + 1)

    @ivy nray_max_backward = maximum(nray[ii, j0, kk])
    @ivy nray_max_forward = maximum(nray[ii, j1, kk])

    nray_max_backward = MPI.Allreduce(nray_max_backward, max, comm)
    nray_max_forward = MPI.Allreduce(nray_max_forward, max, comm)

    send_forward = zeros(length(fields), nray_max_forward, nx + 2, nz + 2)
    send_backward = zeros(length(fields), nray_max_backward, nx + 2, nz + 2)

    receive_backward = zeros(length(fields), nray_max_forward, nx + 2, nz + 2)
    receive_forward = zeros(length(fields), nray_max_backward, nx + 2, nz + 2)

    @ivy for (index, field) in enumerate(fields)
        @. send_forward[index, :, :, :] =
            $getfield(rays, field)[1:nray_max_forward, ii, j1, kk]
        @. send_backward[index, :, :, :] =
            $getfield(rays, field)[1:nray_max_backward, ii, j0, kk]
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

    @ivy for (index, field) in enumerate(fields)
        @. $getfield(rays, field)[1:nray_max_forward, ii, j0 - 1, kk] =
            receive_backward[index, :, :, :]
        @. $getfield(rays, field)[1:nray_max_backward, ii, j1 + 1, kk] =
            receive_forward[index, :, :, :]
    end

    return
end
