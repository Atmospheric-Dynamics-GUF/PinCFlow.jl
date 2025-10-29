"""
```julia
TracerWKBImpact{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
}
```

Container for the gravity-wave-induced tracer fluxes and resulting tracer tendency.

```julia
TracerWKBImpact(nxi::Integer, nyi::Integer, nzi::Integer)::TracerWKBImpact
```

Construct a `TracerWKBImpact` instance with array dimensions given by `nxi`, `nyi`, and `nzi`.

# Fields

  - `uchi::A`: Zonal tracer fluxes due to unresolved gravity waves.

  - `vchi::B`: Meridional tracer fluxes due to unresolved gravity waves.

  - `wchi::C`: Vertical tracer fluxes due to unresolved gravity waves.

  - `dchidt::D`: Leading-order tracer impact of unresolved gravity waves.

# Arguments:

  - `nxi`: Grid-points in `\\widehat{x}`-direction.

  - `nyi`: Grid-points in `\\widehat{y}`-direction.

  - `nzi`: Grid-points in `\\widehat{z}`-direction.
"""
struct TracerWKBImpact{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
}
    uchi::A
    vchi::B
    wchi::C
    dchidt::D
end

function TracerWKBImpact(
    nxi::Integer,
    nyi::Integer,
    nzi::Integer,
)::TracerWKBImpact
    return TracerWKBImpact([zeros(nxi, nyi, nzi) for i in 1:4]...)
end
