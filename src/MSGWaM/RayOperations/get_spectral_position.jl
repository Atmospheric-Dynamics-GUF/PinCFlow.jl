"""
```julia
get_spectral_position(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
```

Return the spectral position of the ray volume specified by ``\\left(r, i, j, k\\right)`` as the tuple `(kr, lr, mr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.
"""
function get_spectral_position end

function get_spectral_position(
    rays::Rays,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{3, <:AbstractFloat}
    @ivy return (rays.k[r, i, j, k], rays.l[r, i, j, k], rays.m[r, i, j, k])
end
