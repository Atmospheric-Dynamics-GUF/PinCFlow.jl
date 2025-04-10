function copy_ray_volume!(
    rays::Rays,
    source::NTuple{4, <:Integer},
    target::NTuple{4, <:Integer},
)
    (irs, ixs, jys, kzs) = source
    (irt, ixt, jyt, kzt) = target

    for field in fieldnames(Rays)
        getfield(rays, field)[irt, ixt, jyt, kzt] =
            getfield(rays, field)[irs, ixs, jys, kzs]
    end

    return
end
