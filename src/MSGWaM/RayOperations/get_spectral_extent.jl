"""
```julia
get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the spectral extents of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.
  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::Float64`: Extent in ``k``-direction.
  - `::Float64`: Extent in ``l``-direction.
  - `::Float64`: Extent in ``m``-direction.
"""
function get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dkray[indices...],
        rays.dlray[indices...],
        rays.dmray[indices...],
    )
end
