"""
```julia
get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the spectral position of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.
  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::Float64`: Position in ``k``-direction.
  - `::Float64`: Position in ``l``-direction.
  - `::Float64`: Position in ``m``-direction.
"""
function get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.k[indices...], rays.l[indices...], rays.m[indices...])
end
