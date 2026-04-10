"""
```julia
Rays{A <: AbstractArray{<:AbstractFloat, 5}}
```

Container for prognostic ray-volume properties.

The field `data` is an array, the dimensions of which represent (in order) the ray-volume properties, the ray volumes within a grid cell, and the three dimensions of the physical grid. The elements of this array can be accessed in two ways:

  - Indexing the field, e.g., `rays.data[1, 1, 1, 1, 1]`.

  - Extracting a property, e.g., `rays.x[1, 1, 1, 1]`.

The property names for the second way are the following:

  - `x`: Position in ``x``.

  - `y`: Position in ``y``.

  - `z`: Position in ``z``.

  - `k`: Position in ``k``.

  - `l`: Position in ``l``.

  - `m`: Position in ``m``.

  - `dxray`: Extent in ``x``.

  - `dyray`: Extent in ``y``.

  - `dzray`: Extent in ``z``.

  - `dkray`: Extent in ``k``.

  - `dlray`: Extent in ``l``.

  - `dmray`: Extent in ``m``.

  - `dens`: Phase-space wave-action density.

```julia
Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)::Rays
```

Construct a `Rays` instance, with an array sized according to the given dimensions.

# Fields

  - `data::A`: Ray-volume data.

# Arguments

  - `nray_wrk`: Size of the spectral dimension of ray-volume arrays.

  - `nxx`: Number of subdomain grid points in ``\\hat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\hat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\hat{z}``-direction.
"""
struct Rays{A <: AbstractArray{<:AbstractFloat, 5}}
    data::A
end

function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)::Rays
    return Rays(zeros(13, nray_wrk, nxx, nyy, nzz))
end
