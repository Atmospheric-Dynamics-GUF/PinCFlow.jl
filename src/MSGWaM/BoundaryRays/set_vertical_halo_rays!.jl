"""
```julia
set_vertical_halo_rays!(state::State)
```

Exchange ray volumes in vertical halo cells.

Performs MPI communication between downward and upward neighbor processes. The number of communicated ray volumes is determined from the maximum counts at the downward and upward boundaries of the MPI subdomains. Solid walls are assumed at the vertical boundaries of the full domain. The corresponding ghost-cell ray volumes are not changed.

# Arguments

  - `state`: Model state.
"""
function set_vertical_halo_rays! end

function set_vertical_halo_rays!(state::State)
    (; comm, zz_size, nzz, nx, ny, ko, i0, i1, j0, j1, k0, k1, down, up) =
        state.domain
    (; nray, rays) = state.wkb

    fields = fieldcount(Rays)

    ii = (i0 - 1):(i1 + 1)
    jj = (j0 - 1):(j1 + 1)

    @ivy nray_max_down = maximum(nray[ii, jj, k0])
    @ivy nray_max_up = maximum(nray[ii, jj, k1])

    nray_max_down = MPI.Allreduce(nray_max_down, max, comm)
    nray_max_up = MPI.Allreduce(nray_max_up, max, comm)

    send_up = zeros(fields, nray_max_up, nx + 2, ny + 2)
    send_down = zeros(fields, nray_max_down, nx + 2, ny + 2)

    receive_down = zeros(fields, nray_max_up, nx + 2, ny + 2)
    receive_up = zeros(fields, nray_max_down, nx + 2, ny + 2)

    @ivy if ko == 0
        for field in 1:fields
            send_up[field, :, :, :] .=
                getfield(rays, field)[1:nray_max_up, ii, jj, k1]
        end

        MPI.Sendrecv!(send_up, receive_up, comm; dest = up, source = up)

        for field in 1:fields
            getfield(rays, field)[1:nray_max_down, ii, jj, k1 + 1] .=
                receive_up[field, :, :, :]
        end
    elseif ko + nzz == zz_size
        for field in 1:fields
            send_down[field, :, :, :] .=
                getfield(rays, field)[1:nray_max_down, ii, jj, k0]
        end

        MPI.Sendrecv!(send_down, receive_down, comm; dest = down, source = down)

        for field in 1:fields
            getfield(rays, field)[1:nray_max_up, ii, jj, k0 - 1] .=
                receive_down[field, :, :, :]
        end
    else
        for field in 1:fields
            send_up[field, :, :, :] .=
                getfield(rays, field)[1:nray_max_up, ii, jj, k1]
            send_down[field, :, :, :] .=
                getfield(rays, field)[1:nray_max_down, ii, jj, k0]
        end

        MPI.Sendrecv!(send_up, receive_down, comm; dest = up, source = down)

        MPI.Sendrecv!(send_down, receive_up, comm; dest = down, source = up)

        for field in 1:fields
            getfield(rays, field)[1:nray_max_up, ii, jj, k0 - 1] .=
                receive_down[field, :, :, :]
            getfield(rays, field)[1:nray_max_down, ii, jj, k1 + 1] .=
                receive_up[field, :, :, :]
        end
    end

    return
end
