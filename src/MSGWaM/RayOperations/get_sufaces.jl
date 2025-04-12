function get_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
    dxr = rays.dxray[indices...]
    dkr = rays.dkray[indices...]
    dyr = rays.dyray[indices...]
    dlr = rays.dlray[indices...]
    dzr = rays.dzray[indices...]
    dmr = rays.dmray[indices...]
    return (dxr * dkr, dyr * dlr, dzr * dmr)
end
