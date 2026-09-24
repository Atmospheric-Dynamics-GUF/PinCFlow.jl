"""
```julia
copy_rays!(
    rays::Rays,
    r::Pair{<:Integer, <:Integer},
    i::Pair{<:Integer, <:Integer},
    j::Pair{<:Integer, <:Integer},
    k::Pair{<:Integer, <:Integer},
)
```

Copy all properties of the ray volume specified by the first components of the index pairs to that specified by the second components.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `r`: Ray-volume indices.

  - `i`: Zonal grid-cell indices.

  - `j`: Meridional grid-cell indices.

  - `k`: Vertical grid-cell indices.
"""
function copy_rays! end

@ivy function copy_rays!(
    rays::Rays,
    r::Pair{<:Integer, <:Integer},
    i::Pair{<:Integer, <:Integer},
    j::Pair{<:Integer, <:Integer},
    k::Pair{<:Integer, <:Integer},
)
    rays.data[:, r[2], i[2], j[2], k[2]] .= rays.data[:, r[1], i[1], j[1], k[1]]

    return
end
