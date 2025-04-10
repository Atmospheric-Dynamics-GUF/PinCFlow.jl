function get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dxray[indices...],
        rays.dyray[indices...],
        rays.dzray[indices...],
    )
end