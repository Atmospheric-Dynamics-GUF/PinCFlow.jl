function get_physical_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (rays.x[indices...], rays.y[indices...], rays.z[indices...])
end