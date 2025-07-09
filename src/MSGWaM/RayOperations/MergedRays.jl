"""
```julia
MergedRays{A <: AbstractVector{<:AbstractFloat}, B <: Ref{<:AbstractFloat}}
```

Composite type used for creating merged ray volumes.

# Fields

  - `xr::A`: Outermost ray-volume bounds in ``x``.
  - `yr::A`: Outermost ray-volume bounds in ``y``.
  - `zr::A`: Outermost ray-volume bounds in ``z``.
  - `kr::A`: Outermost ray-volume bounds in ``k``.
  - `lr::A`: Outermost ray-volume bounds in ``l``.
  - `mr::A`: Outermost ray-volume bounds in ``m``.
  - `nr::B`: Wave-action integral.
"""
struct MergedRays{
    A <: AbstractVector{<:AbstractFloat},
    B <: Ref{<:AbstractFloat},
}
    xr::A
    yr::A
    zr::A
    kr::A
    lr::A
    mr::A
    nr::B
end
