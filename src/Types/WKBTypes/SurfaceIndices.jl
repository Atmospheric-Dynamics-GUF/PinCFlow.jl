"""
```julia
SurfaceIndices{A <: AbstractArray{<:Integer, 3}, B <: AbstractVector{<:Integer}}
```

Indices that connect orographic wave modes to the corresponding ray volumes launched by [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref).

```julia
SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)
```

Construct a `SurfaceIndices` instance, with arrays sized according to the given dimensions.

# Fields

- `ir_sfc::A`: Ray-volume indices.
- `ix2_sfc::B`: Zonal indices within grid cells.
- `jy2_sfc::B`: Meridional indices within grid cells.
- `kz2_sfc::B`: Vertical indices within grid cells.
- `ik_sfc::B`: Index in ``k``.
- `jl_sfc::B`: Index in ``l``.
- `km_sfc::B`: Index in ``m``.
- `iwm_sfc::B`: Wave-mode index.

# Arguments

- `n_sfc`: Number of orographic wave modes per grid cell.
- `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.
- `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
"""
struct SurfaceIndices{
    A <: AbstractArray{<:Integer, 3},
    B <: AbstractVector{<:Integer},
}
    ir_sfc::A
    ix2_sfc::B
    jy2_sfc::B
    kz2_sfc::B
    ik_sfc::B
    jl_sfc::B
    km_sfc::B
    iwm_sfc::B
end

function SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)
    return SurfaceIndices(
        zeros(Int, n_sfc, nxx, nyy),
        [zeros(Int, n_sfc) for i in 1:7]...,
    )
end
