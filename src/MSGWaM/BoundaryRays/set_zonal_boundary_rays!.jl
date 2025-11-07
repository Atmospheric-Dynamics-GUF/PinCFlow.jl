"""
```julia
set_zonal_boundary_rays!(state::State)
```

Set the ray volumes at the zonal boundaries.

This method first enforces zonal boundary conditions for `state.wkb.nray` (by applying `set_zonal_boundaries_of_field!` to it) and then sets the corresponding boundary ray volumes, assuming periodicity. If the domain is parallelized in ``\\widehat{x}``, ray volumes are communicated between MPI processes, using `set_zonal_halo_rays!`. At the zonal boundaries of the domain, the ``x``-coordinates of ray volumes are adjusted such that shifting works properly.

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_zonal_halo_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function set_zonal_boundary_rays! end

function set_zonal_boundary_rays!(state::State)
    (; namelists, domain) = state
    (; x_size, npx) = namelists.domain
    (; nx, io, i0, i1, j0, j1, k0, k1) = domain
    (; lx, x) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count.
    set_zonal_boundaries_of_field!(nray, namelists, domain; layers = (1, 1, 1))

    # Set ray-volumes properties.
    @ivy if npx > 1
        set_zonal_halo_rays!(state)
    else
        for k in (k0 - 1):(k1 + 1), j in (j0 - 1):(j1 + 1)
            for r in 1:nray[i0 - 1, j, k]
                copy_rays!(rays, r => r, i1 => i0 - 1, j => j, k => k)
            end

            for r in 1:nray[i1 + 1, j, k]
                copy_rays!(rays, r => r, i0 => i1 + 1, j => j, k => k)
            end
        end
    end

    @ivy if io == 0
        for k in (k0 - 1):(k1 + 1), j in (j0 - 1):(j1 + 1), i in (i0 - 1):i0
            for r in 1:nray[i, j, k]
                xr = rays.x[r, i, j, k]
                xrt = xr - lx

                if abs(xrt - x[i]) < abs(xr - x[i])
                    xr = xrt
                end

                rays.x[r, i, j, k] = xr
            end
        end
    end

    @ivy if io + nx == x_size
        for k in (k0 - 1):(k1 + 1), j in (j0 - 1):(j1 + 1), i in i1:(i1 + 1)
            for r in 1:nray[i, j, k]
                xr = rays.x[r, i, j, k]
                xrt = xr + lx

                if abs(xrt - x[i]) < abs(xr - x[i])
                    xr = xrt
                end

                rays.x[r, i, j, k] = xr
            end
        end
    end

    return
end
