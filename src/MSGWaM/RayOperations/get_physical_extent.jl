"""
```julia
get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the physical extents of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.
  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::Float64`: Extent in ``x``-direction.
  - `::Float64`: Extent in ``y``-direction.
  - `::Float64`: Extent in ``z``-direction.
"""
function get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dxray[indices...],
        rays.dyray[indices...],
        rays.dzray[indices...],
    )
end
