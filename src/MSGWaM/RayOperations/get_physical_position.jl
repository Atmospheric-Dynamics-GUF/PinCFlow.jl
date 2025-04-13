function get_physical_position(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        copy(rays.x[indices...]),
        copy(rays.y[indices...]),
        copy(rays.z[indices...]),
    )
end
