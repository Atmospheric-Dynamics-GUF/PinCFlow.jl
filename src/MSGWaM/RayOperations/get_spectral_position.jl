function get_spectral_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        copy(rays.k[indices...]),
        copy(rays.l[indices...]),
        copy(rays.m[indices...]),
    )
end
