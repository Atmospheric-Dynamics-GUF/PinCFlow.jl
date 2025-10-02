"""
```julia
set_vertical_boundary_rays!(state::State)
```

Enforce vertical boundary conditions for ray volumes.

If the domain is parallelized in ``\\widehat{z}``, ray-volume counts and the ray volumes themselves are first communicated between MPI processes, using `set_vertical_halos_of_field!` and `set_vertical_halo_rays!`, respectively. The vertical boundary conditions are then enforced by cutting (removing) ray volumes that have partially (fully) crossed the upper boundary and reflecting ray volumes (by adjusting the vertical position and wavenumber) that have at least partially crossed the lower boundary from above.

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.MPIOperations.set_vertical_halos_of_field!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_vertical_halo_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function set_vertical_boundary_rays! end

function set_vertical_boundary_rays!(state::State)
    (; namelists, domain) = state
    (; npz) = namelists.domain
    (; zz_size, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; lx, ly, lz, dx, dy, topography_surface) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count and ray-volumes properties.
    if npz > 1
        set_vertical_halos_of_field!(
            nray,
            namelists,
            domain;
            layers = (1, 1, 1),
        )
        set_vertical_halo_rays!(state)
    end

    # Reflect ray volumes at the lower boundary.
    @ivy if ko == 0
        for k in k0:(k0 + 1), j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
            for r in 1:nray[i, j, k]
                xr = rays.x[r, i, j, k]
                yr = rays.y[r, i, j, k]
                zr = rays.z[r, i, j, k]
                dzr = rays.dzray[r, i, j, k]
                wnrm = rays.m[r, i, j, k]

                iray = floor(Int, (xr + lx / 2) / dx) + i0 - io
                jray = floor(Int, (yr + ly / 2) / dy) + j0 - jo
                if topography_surface[iray, jray] - zr + 0.5 * dzr > eps()
                    rays.z[r, i, j, k] =
                        2.0 * topography_surface[iray, jray] - zr + dzr
                    rays.m[r, i, j, k] = -wnrm
                end
            end
        end
    end

    # Cut ray volumes at the upper boundary.
    @ivy if ko + nzz == zz_size
        for k in (k1 - 1):k1, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
            local_count = 0
            for r in 1:nray[i, j, k]
                zr = rays.z[r, i, j, k]
                dzr = rays.dzray[r, i, j, k]

                if zr - 0.5 * dzr > lz
                    continue
                end
                if zr + 0.5 * dzr > lz
                    rays.dzray[r, i, j, k] = lz - zr + 0.5 * dzr
                    rays.z[r, i, j, k] = lz - 0.5 * rays.dzray[r, i, j, k]
                end

                local_count += 1
                if local_count != r
                    copy_rays!(rays, r => local_count, i => i, j => j, k => k)
                    rays.dens[r, i, j, k] = 0.0
                end
            end
            nray[i, j, k] = local_count
        end
    end

    return
end
