function get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.k[indices...], rays.l[indices...], rays.m[indices...])
end