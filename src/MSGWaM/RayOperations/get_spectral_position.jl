"""
```julia
get_spectral_position(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
```

Return the spectral position of the ray volume specified by `indices` as the tuple `(kr, lr, mr)`.

# Arguments

  - `rays`: Collection of ray-volume-property arrays.

  - `indices`: Indices of the ray volume of interest.
"""
function get_spectral_position end

function get_spectral_position(
    rays::Rays,
    indices::NTuple{4, <:Integer},
)::NTuple{3, <:AbstractFloat}
    @ivy return (rays.k[indices...], rays.l[indices...], rays.m[indices...])
end
