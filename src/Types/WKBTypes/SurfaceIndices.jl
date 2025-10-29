"""
```julia
SurfaceIndices{
    A <: AbstractArray{<:Integer, 3},
    B <: AbstractVector{<:Integer},
    C <: AbstractVector{<:Integer},
    D <: AbstractVector{<:Integer},
    E <: AbstractVector{<:Integer},
    F <: AbstractVector{<:Integer},
    G <: AbstractVector{<:Integer},
    H <: AbstractVector{<:Integer},
}
```

Indices that connect orographic wave modes to the corresponding ray volumes launched by [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref).

```julia
SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)::SurfaceIndices
```

Construct a `SurfaceIndices` instance, with arrays sized according to the given dimensions.

# Fields

  - `rs::A`: Ray-volume indices.

  - `ixs::B`: Zonal indices within grid cells.

  - `jys::C`: Meridional indices within grid cells.

  - `kzs::D`: Vertical indices within grid cells.

  - `iks::E`: Indices in ``k``.

  - `jls::F`: Indices in ``l``.

  - `kms::G`: Indices in ``m``.

  - `alphas::H`: Wave-mode indices.

# Arguments

  - `n_sfc`: Number of orographic wave modes per grid cell.

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
"""
struct SurfaceIndices{
    A <: AbstractArray{<:Integer, 3},
    B <: AbstractVector{<:Integer},
    C <: AbstractVector{<:Integer},
    D <: AbstractVector{<:Integer},
    E <: AbstractVector{<:Integer},
    F <: AbstractVector{<:Integer},
    G <: AbstractVector{<:Integer},
    H <: AbstractVector{<:Integer},
}
    rs::A
    ixs::B
    jys::C
    kzs::D
    iks::E
    jls::F
    kms::G
    alphas::H
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
