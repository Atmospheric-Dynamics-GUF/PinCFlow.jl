"""
```julia
get_spectral_extent(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
```

Return the spectral extents of the ray volume specified by ``\\left(r, i, j, k\\right)`` as the tuple `(dkr, dlr, dmr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.
"""
function get_spectral_extent end

function get_spectral_extent(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
    @ivy return (
        rays.dkray[r, i, j, k],
        rays.dlray[r, i, j, k],
        rays.dmray[r, i, j, k],
    )
end
