"""
```julia
get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})::NTuple{3, <:AbstractFloat}
```

Return the physical extents of the ray volume specified by `indices` as tuple `(dx, dy, dz)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `indices`: Indices of the ray volume of interest.
"""
function get_physical_extent end

function get_physical_extent(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
    return (
        rays.dxray[indices...],
        rays.dyray[indices...],
        rays.dzray[indices...],
    )
end
