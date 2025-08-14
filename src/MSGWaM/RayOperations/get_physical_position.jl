"""
```julia
get_physical_position(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
```

Return the physical position of the ray volume specified by `indices` as the tuple `(xr, yr, zr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `indices`: Indices of the ray volume of interest.
"""
function get_physical_position end

function get_physical_position(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
    return (rays.x[indices...], rays.y[indices...], rays.z[indices...])
end
