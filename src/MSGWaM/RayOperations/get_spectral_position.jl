"""
```julia
get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
```

Return the spectral position of the ray volume specified by `indices`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.
  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::AbstractFloat`: Position in ``k``-direction.
  - `::AbstractFloat`: Position in ``l``-direction.
  - `::AbstractFloat`: Position in ``m``-direction.
"""
function get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.k[indices...], rays.l[indices...], rays.m[indices...])
end
