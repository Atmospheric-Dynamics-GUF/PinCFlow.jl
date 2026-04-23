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
    (; z_size) = state.namelists.domain
    (; comm, nz, nx, ny, ko, i0, i1, j0, j1, k0, k1, down, up) = state.domain
    (; nray, rays) = state.wkb

    ii = (i0 - 1):(i1 + 1)
    jj = (j0 - 1):(j1 + 1)

    @ivy nray_max_down = maximum(nray[ii, jj, k0])
    @ivy nray_max_up = maximum(nray[ii, jj, k1])

    nray_max_down = MPI.Allreduce(nray_max_down, max, comm)
    nray_max_up = MPI.Allreduce(nray_max_up, max, comm)

    @ivy if ko == 0
        if nray_max_up > 0
            MPI.Send(rays.data[:, 1:nray_max_up, ii, jj, k1], comm; dest = up)
        end

        if nray_max_down > 0
            MPI.Recv!(
                rays.data[:, 1:nray_max_down, ii, jj, k1 + 1],
                comm;
                source = up,
            )
        end
    elseif ko + nz == z_size
        if nray_max_up > 0
            MPI.Recv!(
                rays.data[:, 1:nray_max_up, ii, jj, k0 - 1],
                comm;
                source = down,
            )
        end

        if nray_max_down > 0
            MPI.Send(
                rays.data[:, 1:nray_max_down, ii, jj, k0],
                comm;
                dest = down,
            )
        end
    else
        if nray_max_up > 0
            MPI.Sendrecv!(
                rays.data[:, 1:nray_max_up, ii, jj, k1],
                rays.data[:, 1:nray_max_up, ii, jj, k0 - 1],
                comm;
                dest = up,
                source = down,
            )
        end

        if nray_max_down > 0
            MPI.Sendrecv!(
                rays.data[:, 1:nray_max_down, ii, jj, k0],
                rays.data[:, 1:nray_max_down, ii, jj, k1 + 1],
                comm;
                dest = down,
                source = up,
            )
        end
    end

    return
end
