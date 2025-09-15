"""
```julia
get_physical_extent(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
```

Return the physical extents of the ray volume specified by ``\\left(r, i, j, k\\right)`` as the tuple `(dxr, dyr, dzr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.
"""
function get_physical_extent end

function get_physical_extent(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
    @ivy return (
        rays.dxray[r, i, j, k],
        rays.dyray[r, i, j, k],
        rays.dzray[r, i, j, k],
    )
end
