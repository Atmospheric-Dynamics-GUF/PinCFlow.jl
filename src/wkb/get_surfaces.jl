function get_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.area_xk[indices...],
        rays.area_yl[indices...],
        rays.area_zm[indices...],
    )
end