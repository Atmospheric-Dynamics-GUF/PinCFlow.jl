"""
```julia
set_meridional_boundary_rays!(state::State)
```

Set the ray volumes at the meridional boundaries.

This method first enforces meridional boundary conditions for `state.wkb.nray` (by applying `set_meridional_boundaries_of_field!` to it) and then sets the corresponding boundary ray volumes, assuming periodicity. If the domain is parallelized in ``\\widehat{y}``, ray volumes are communicated between MPI processes, using `set_meridional_halo_rays!`. At the meridional boundaries of the domain, the ``y``-coordinates of ray volumes are adjusted such that shifting works properly.

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)

  - [`PinCFlow.MSGWaM.BoundaryRays.set_meridional_halo_rays!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function set_meridional_boundary_rays! end

function set_meridional_boundary_rays!(state::State)
    (; namelists, domain) = state
    (; npy) = namelists.domain
    (; nyy, yy_size, jo, i0, i1, j0, j1, k0, k1) = domain
    (; ly, y) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count.
    set_meridional_boundaries_of_field!(
        nray,
        namelists,
        domain;
        layers = (1, 1, 1),
    )

    # Set ray-volumes properties.
    @ivy if npy > 1
        set_meridional_halo_rays!(state)
    else
        for k in (k0 - 1):(k1 + 1), i in (i0 - 1):(i1 + 1)
            for r in 1:nray[i, j0 - 1, k]
                copy_rays!(rays, r => r, i => i, j1 => j0 - 1, k => k)
            end

            for r in 1:nray[i, j1 + 1, k]
                copy_rays!(rays, r => r, i => i, j0 => j1 + 1, k => k)
            end
        end
    end

    @ivy if jo == 0
        for k in (k0 - 1):(k1 + 1), j in (j0 - 1):j0, i in (i0 - 1):(i1 + 1)
            for r in 1:nray[i, j, k]
                yr = rays.y[r, i, j, k]
                yrt = yr - ly

                if abs(yrt - y[j]) < abs(yr - y[j])
                    yr = yrt
                end

                rays.y[r, i, j, k] = yr
            end
        end
    end

    @ivy if jo + nyy == yy_size
        for k in (k0 - 1):(k1 + 1), j in j1:(j1 + 1), i in (i0 - 1):(i1 + 1)
            for r in 1:nray[i, j, k]
                yr = rays.y[r, i, j, k]
                yrt = yr + ly

                if abs(yrt - y[j]) < abs(yr - y[j])
                    yr = yrt
                end

                rays.y[r, i, j, k] = yr
            end
        end
    end

    return
end
