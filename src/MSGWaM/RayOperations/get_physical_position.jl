"""
```julia
get_physical_position(rays::Rays, indices::NTuple{4, <:Integer})
```

# Arguments

  - `rays::Rays`: Ray volume data structure
  - `indices::NTuple{4, <:Integer}`: Ray indices (iray, ix, jy, kz)
"""
function get_physical_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.x[indices...], rays.y[indices...], rays.z[indices...])
end
