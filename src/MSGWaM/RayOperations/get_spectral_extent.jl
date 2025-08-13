"""
```julia
get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the spectral extents of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::AbstractFloat`: Extent in ``k``-direction.

  - `::AbstractFloat`: Extent in ``l``-direction.

  - `::AbstractFloat`: Extent in ``m``-direction.
"""
function get_spectral_extent end

function get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dkray[indices...],
        rays.dlray[indices...],
        rays.dmray[indices...],
    )
end
