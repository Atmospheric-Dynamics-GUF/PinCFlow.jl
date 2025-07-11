"""
```julia
get_physical_position(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the physical position of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.
  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::Float64`: Position in ``x``-direction.
  - `::Float64`: Position in ``y``-direction.
  - `::Float64`: Position in ``z``-direction.
"""
function get_physical_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.x[indices...], rays.y[indices...], rays.z[indices...])
end
