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
function set_zonal_boundary_rays!(state::State)
    (; namelists, domain) = state
    (; npx) = namelists.domain
    (; nxx, sizexx, io, i0, i1, j0, j1, k0, k1) = domain
    (; lx, x) = state.grid
    (; nray, rays) = state.wkb

    # Set ray-volume count.
    set_zonal_boundaries_of_field!(nray, namelists, domain; layers = (1, 1, 1))

    # Set ray-volumes properties.
    if npx > 1
        set_zonal_halo_rays!(state)
    else
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1)
            for iray in 1:nray[i0 - 1, jy, kz]
                copy_rays!(rays, (iray, i1, jy, kz), (iray, i0 - 1, jy, kz))
            end

            for iray in 1:nray[i1 + 1, jy, kz]
                copy_rays!(rays, (iray, i0, jy, kz), (iray, i1 + 1, jy, kz))
            end
        end
    end

    if io == 0
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):i0
            for iray in 1:nray[ix, jy, kz]
                xr = rays.x[iray, ix, jy, kz]
                xrt = xr - lx[2] + lx[1]

                if abs(xrt - x[io + ix]) < abs(xr - x[io + ix])
                    xr = xrt
                end

                rays.x[iray, ix, jy, kz] = xr
            end
        end
    end

    if io + nxx == sizexx
        for kz in (k0 - 1):(k1 + 1), jy in (j0 - 1):(j1 + 1), ix in i1:(i1 + 1)
            for iray in 1:nray[ix, jy, kz]
                xr = rays.x[iray, ix, jy, kz]
                xrt = xr + lx[2] - lx[1]

                if abs(xrt - x[io + ix]) < abs(xr - x[io + ix])
                    xr = xrt
                end

                rays.x[iray, ix, jy, kz] = xr
            end
        end
    end

    return
end
