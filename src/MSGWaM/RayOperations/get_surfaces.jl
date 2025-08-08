"""
```julia
get_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
```

Compute phase-space surface areas of the ray volume specified by `indices`.

# Arguments

- `rays`: Collection of ray-volume-property arrays.
- `indices`: Indices of the ray volume of interest.

# Returns

- `::AbstractFloat`: Surface area in ``x``-``k`` subspace.
- `::AbstractFloat`: Surface area in ``y``-``l`` subspace.
- `::AbstractFloat`: Surface area in ``z``-``m`` subspace.
"""
function get_surfaces end

function get_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
    (dxr, dyr, dzr) = get_physical_extent(rays, indices)
    (dkr, dlr, dmr) = get_spectral_extent(rays, indices)
    return (dxr * dkr, dyr * dlr, dzr * dmr)
end
