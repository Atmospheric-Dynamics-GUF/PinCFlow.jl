"""
```julia
remove_rays!(state::State)
```

Remove gaps (i.e. zero-wave-action ray volumes between nonzero-wave-action ray volumes) in the ray-volume arrays.

In each grid cell, this method moves all ray volumes as far to the front of the arrays possible and updates `nray` accordingly, so that every ray volume in the range `1:nray[i, j, k]` has nonzero wave action.

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function remove_rays! end

function remove_rays!(state::State)
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray, rays) = state.wkb

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == sizezz ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
        local_count = 0
        for r in 1:nray[i, j, k]
            if rays.dens[r, i, j, k] == 0
                continue
            end
            local_count += 1
            if local_count != r
                copy_rays!(rays, r => local_count, i => i, j => j, k => k)
                rays.dens[r, i, j, k] = 0.0
            end
        end
        nray[i, j, k] = local_count
    end

    return
end
