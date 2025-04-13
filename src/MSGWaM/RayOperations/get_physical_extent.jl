function get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        copy(rays.dxray[indices...]),
        copy(rays.dyray[indices...]),
        copy(rays.dzray[indices...]),
    )
end
