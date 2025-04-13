function get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        copy(rays.dkray[indices...]),
        copy(rays.dlray[indices...]),
        copy(rays.dmray[indices...]),
    )
end
