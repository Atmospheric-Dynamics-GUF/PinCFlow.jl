"""
```julia
get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the physical extents of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.
  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::AbstractFloat`: Extent in ``x``-direction.
  - `::AbstractFloat`: Extent in ``y``-direction.
  - `::AbstractFloat`: Extent in ``z``-direction.
"""
function get_physical_extent end

function get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dxray[indices...],
        rays.dyray[indices...],
        rays.dzray[indices...],
    )
end
