"""
```julia
MergedRays{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractVector{<:AbstractFloat},
}
```

Composite type used for creating merged ray volumes.

```julia
MergedRays(bounds::Integer, count::Integer)::MergedRays
```

Construct a `MergedRays` instance, with arrays sized according to the given dimensions.

# Fields

  - `xr::A`: Outermost ray-volume bounds in ``x``.

  - `yr::A`: Outermost ray-volume bounds in ``y``.

  - `zr::A`: Outermost ray-volume bounds in ``z``.

  - `kr::A`: Outermost ray-volume bounds in ``k``.

  - `lr::A`: Outermost ray-volume bounds in ``l``.

  - `mr::A`: Outermost ray-volume bounds in ``m``.

  - `nr::B`: Wave-action integral.

# Arguments

  - `bounds`: Number of bounds in each dimension.

  - `count`: Maximum ray-volume count per grid cell.
"""
struct MergedRays{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractVector{<:AbstractFloat},
}
    xr::A
    yr::A
    zr::A
    kr::A
    lr::A
    mr::A
    nr::B
end

function MergedRays(bounds::Integer, count::Integer)::MergedRays
    return MergedRays([zeros(bounds, count) for i in 1:6]..., zeros(count))
end
