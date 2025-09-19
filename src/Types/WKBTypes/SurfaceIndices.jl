"""
```julia
SurfaceIndices{A <: AbstractArray{<:Integer, 3}, B <: AbstractVector{<:Integer}}
```

Indices that connect orographic wave modes to the corresponding ray volumes launched by [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref).

```julia
SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)::SurfaceIndices
```

Construct a `SurfaceIndices` instance, with arrays sized according to the given dimensions.

# Fields

  - `rs::A`: Ray-volume indices.

  - `ixs::B`: Zonal indices within grid cells.

  - `jys::B`: Meridional indices within grid cells.

  - `kzs::B`: Vertical indices within grid cells.

  - `iks::B`: Indices in ``k``.

  - `jls::B`: Indices in ``l``.

  - `kms::B`: Indices in ``m``.

  - `alphas::B`: Wave-mode indices.

# Arguments

  - `n_sfc`: Number of orographic wave modes per grid cell.

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
"""
struct SurfaceIndices{
    A <: AbstractArray{<:Integer, 3},
    B <: AbstractVector{<:Integer},
}
    rs::A
    ixs::B
    jys::B
    kzs::B
    iks::B
    jls::B
    kms::B
    alphas::B
end

function SurfaceIndices(
    n_sfc::Integer,
    nxx::Integer,
    nyy::Integer,
)::SurfaceIndices
    return SurfaceIndices(
        zeros(Int, n_sfc, nxx, nyy),
        [zeros(Int, n_sfc) for i in 1:7]...,
    )
end
