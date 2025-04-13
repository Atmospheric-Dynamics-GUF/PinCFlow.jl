function get_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
    (dxr, dyr, dzr) = get_physical_extent(rays, indices)
    (dkr, dlr, dmr) = get_spectral_extent(rays, indices)
    return (dxr * dkr, dyr * dlr, dzr * dmr)
end
