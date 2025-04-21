function get_spectral_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dkray[indices...],
        rays.dlray[indices...],
        rays.dmray[indices...],
    )
end
