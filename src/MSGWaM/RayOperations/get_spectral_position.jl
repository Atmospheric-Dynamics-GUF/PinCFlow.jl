"""
    get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})

# Arguments

  - `rays::Rays`: Ray volume data structure
  - `indices::NTuple{4, <:Integer}`: Ray indices (iray, ix, jy, kz)
"""
function get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.k[indices...], rays.l[indices...], rays.m[indices...])
end
