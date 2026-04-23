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

    jj = (j0 - 1):(j1 + 1)
    kk = (k0 - 1):(k1 + 1)

    @ivy nray_max_left = maximum(nray[i0, jj, kk])
    @ivy nray_max_right = maximum(nray[i1, jj, kk])

    nray_max_left = MPI.Allreduce(nray_max_left, max, comm)
    nray_max_right = MPI.Allreduce(nray_max_right, max, comm)

    @ivy if nray_max_right > 0
        MPI.Sendrecv!(
            rays.data[:, 1:nray_max_right, i1, jj, kk],
            rays.data[:, 1:nray_max_right, i0 - 1, jj, kk],
            comm;
            dest = right,
            source = left,
        )
    end

    @ivy if nray_max_left > 0
        MPI.Sendrecv!(
            rays.data[:, 1:nray_max_left, i0, jj, kk],
            rays.data[:, 1:nray_max_left, i1 + 1, jj, kk],
            comm;
            dest = left,
            source = right,
        )
    end

    return
end
