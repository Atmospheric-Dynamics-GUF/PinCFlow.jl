"""
```julia
get_surfaces(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
```

Compute phase-space surface areas of the ray volume specified by ``\\left(r, i, j, k\\right)`` and return them as the tuple `(axk, ayl, azm)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.
"""
function get_surfaces end

function get_surfaces(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
    (dxr, dyr, dzr) = get_physical_extent(rays, r, i, j, k)
    (dkr, dlr, dmr) = get_spectral_extent(rays, r, i, j, k)
    return (dxr * dkr, dyr * dlr, dzr * dmr)
end
