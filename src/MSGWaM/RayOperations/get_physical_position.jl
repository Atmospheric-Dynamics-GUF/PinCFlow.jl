"""
```julia
get_physical_position(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
```

Return the physical position of the ray volume specified by ``\\left(r, i, j, k\\right)`` as the tuple `(xr, yr, zr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.
"""
function get_physical_position end

function get_physical_position(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
    @ivy return (rays.x[r, i, j, k], rays.y[r, i, j, k], rays.z[r, i, j, k])
end
