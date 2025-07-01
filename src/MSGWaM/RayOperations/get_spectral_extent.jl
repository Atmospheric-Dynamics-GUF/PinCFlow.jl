"""
```julia
get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
```

# Arguments

  - `rays::Rays`: Ray volume data structure
  - `indices::NTuple{4, <:Integer}`: Ray indices (iray, ix, jy, kz)
"""
function get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dkray[indices...],
        rays.dlray[indices...],
        rays.dmray[indices...],
    )
end
