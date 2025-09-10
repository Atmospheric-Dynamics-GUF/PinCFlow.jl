"""
```julia
get_spectral_extent(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
```

Return the spectral extents of the ray volume specified by `indices` as the tuple `(dkr, dlr, dmr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `indices`: Indices of the ray volume of interest.
"""
function get_spectral_extent end

function get_spectral_extent(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
    @ivy return (
        rays.dkray[indices...],
        rays.dlray[indices...],
        rays.dmray[indices...],
    )
end
