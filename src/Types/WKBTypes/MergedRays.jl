"""
```julia
MergedRays{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractMatrix{<:AbstractFloat},
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractMatrix{<:AbstractFloat},
    E <: AbstractMatrix{<:AbstractFloat},
    F <: AbstractMatrix{<:AbstractFloat},
    G <: AbstractVector{<:AbstractFloat},
}
```

Composite type used for creating merged ray volumes.

```julia
MergedRays(bounds::Integer, count::Integer)::MergedRays
```

Construct a `MergedRays` instance, with arrays sized according to the given dimensions.

# Fields

  - `xr::A`: Outermost ray-volume bounds in ``x``.

  - `yr::B`: Outermost ray-volume bounds in ``y``.

  - `zr::C`: Outermost ray-volume bounds in ``z``.

  - `kr::D`: Outermost ray-volume bounds in ``k``.

  - `lr::E`: Outermost ray-volume bounds in ``l``.

  - `mr::F`: Outermost ray-volume bounds in ``m``.

  - `nr::G`: Wave-action integral.

# Arguments

  - `bounds`: Number of bounds in each dimension.

  - `count`: Maximum ray-volume count per grid cell.
"""
struct MergedRays{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractMatrix{<:AbstractFloat},
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractMatrix{<:AbstractFloat},
    E <: AbstractMatrix{<:AbstractFloat},
    F <: AbstractMatrix{<:AbstractFloat},
    G <: AbstractVector{<:AbstractFloat},
}
    xr::A
    yr::B
    zr::C
    kr::D
    lr::E
    mr::F
    nr::G
end

function MergedRays(bounds::Integer, count::Integer)::MergedRays
    return MergedRays([zeros(bounds, count) for i in 1:6]..., zeros(count))
end
