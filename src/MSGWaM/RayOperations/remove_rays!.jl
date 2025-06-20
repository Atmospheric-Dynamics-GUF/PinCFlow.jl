"""
    remove_rays!(state::State)

Remove ray volumes with zero wave action density and compact ray arrays.

Performs cleanup by removing rays marked for deletion (with zero density) and
compacting the remaining rays to eliminate gaps in the ray volume arrays.
This maintains array efficiency and prevents memory fragmentation.

# Arguments

  - `state::State`: Complete simulation state containing ray data

# Algorithm

 1. **Extended Domain**: Processes rays in extended domain including halo regions
 2. **Compaction**: For each grid cell, compacts non-zero rays to beginning of array
 3. **Array Cleanup**: Moves valid rays to fill gaps left by removed rays
 4. **Count Update**: Updates `nray` count to reflect actual number of active rays

# Grid Cell Range

Processes rays in extended domain:

  - X: `(i0-1):(i1+1)` or adjusted for domain boundaries
  - Y: `(j0-1):(j1+1)`
  - Z: Adjusted based on vertical position in global domain

# Implementation Details

  - Uses `copy_rays!` to move ray data between array positions
  - Sets density to zero for rays in their original positions after copying
  - Maintains ray order while eliminating gaps
  - Efficient O(n) compaction algorithm per grid cell

# Memory Management

Essential for:

  - Preventing memory waste from "dead" ray volumes
  - Maintaining array access efficiency
  - Enabling accurate ray counting for diagnostics
  - Preparing arrays for subsequent operations

# Usage Context

Typically called after:

  - Ray shifting operations
  - Boundary condition applications
  - Saturation scheme applications
  - Any operation that may mark rays for removal
"""
function remove_rays!(state::State)
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; nray, rays) = state.wkb

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for kz in kz0:kz1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        nrlc = 0
        for iray in 1:nray[ix, jy, kz]
            if rays.dens[iray, ix, jy, kz] == 0
                continue
            end
            nrlc += 1
            if nrlc != iray
                copy_rays!(rays, (iray, ix, jy, kz), (nrlc, ix, jy, kz))
                rays.dens[iray, ix, jy, kz] = 0.0
            end
        end
        nray[ix, jy, kz] = nrlc
    end

    return
end
