"""
```julia
copy_rays!(
    rays::Rays,
    source::NTuple{4, <:Integer},
    target::NTuple{4, <:Integer},
)
```

Copy all properties of the ray volume specified by `source` to that specified by `target`.

# Arguments

- `rays`: Collection of ray-volume-property arrays.
- `source`: Indices of the source ray volume.
- `target`: Indices of the target ray volume.
"""
function copy_rays! end

function copy_rays!(
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
