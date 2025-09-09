"""
```julia
set_zonal_halo_rays!(state::State)
```

Exchange ray volumes in zonal halo cells.

Performs bidirectional MPI communication between left and right neighbor processes. The number of communicated ray volumes is determined from the maximum counts at the left and right boundaries of the MPI subdomains.

# Arguments

  - `state`: Model state.
"""
function set_zonal_halo_rays! end

function set_zonal_halo_rays!(state::State)
    (; comm, ny, nz, i0, i1, j0, j1, k0, k1, left, right) = state.domain
    (; nray, rays) = state.wkb

    fields = fieldnames(Rays)

    jj = (j0 - 1):(j1 + 1)
    kk = (k0 - 1):(k1 + 1)

    @views nray_max_left = maximum(nray[i0, jj, kk])
    @views nray_max_right = maximum(nray[i1, jj, kk])

    nray_max_left = MPI.Allreduce(nray_max_left, max, comm)
    nray_max_right = MPI.Allreduce(nray_max_right, max, comm)

    send_right = zeros(length(fields), nray_max_right, ny + 2, nz + 2)
    send_left = zeros(length(fields), nray_max_left, ny + 2, nz + 2)

    receive_left = zeros(length(fields), nray_max_right, ny + 2, nz + 2)
    receive_right = zeros(length(fields), nray_max_left, ny + 2, nz + 2)

    for (index, field) in enumerate(fields)
        @views send_right[index, :, :, :] .=
            getfield(rays, field)[1:nray_max_right, i1, jj, kk]
        @views send_left[index, :, :, :] .=
            getfield(rays, field)[1:nray_max_left, i0, jj, kk]
    end

    MPI.Sendrecv!(send_right, receive_left, comm; dest = right, source = left)

    MPI.Sendrecv!(send_left, receive_right, comm; dest = left, source = right)

    for (index, field) in enumerate(fields)
        @views getfield(rays, field)[1:nray_max_right, i0 - 1, jj, kk] .=
            receive_left[index, :, :, :]
        @views getfield(rays, field)[1:nray_max_left, i1 + 1, jj, kk] .=
            receive_right[index, :, :, :]
    end

    return
end
