"""
```julia
get_surfaces(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
```

Compute phase-space surface areas of the ray volume specified by `indices` and return them as the tuple `(axk, ayl, azm)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `indices`: Indices of the ray volume of interest.
"""
function get_surfaces end

function get_surfaces(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
    (dxr, dyr, dzr) = get_physical_extent(rays, indices)
    (dkr, dlr, dmr) = get_spectral_extent(rays, indices)
    return (dxr * dkr, dyr * dlr, dzr * dmr)
end
