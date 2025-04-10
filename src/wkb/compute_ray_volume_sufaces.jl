function compute_ray_volume_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
    dx = rays.dxray[indices...]
    dk = rays.dkray[indices...]
    dy = rays.dyray[indices...]
    dl = rays.dlray[indices...]
    dz = rays.dzray[indices...]
    dm = rays.dmray[indices...]

    return (dx * dk, dy * dl, dz * dm)
end
