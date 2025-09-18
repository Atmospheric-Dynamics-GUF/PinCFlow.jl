"""
```julia
TracerWKBImpact{A <: AbstractArray{<:AbstractFloat, 3}}
```

Container for the gravity-wave-induced tracer fluxes and resulting tracer tendency.

```julia
TracerWKBImpact(nxi::Integer, nyi::Integer, nzi::Integer)
```

Construct a `TracerWKBImpact` instance with array dimensions given by `nxi`, `nyi`, and `nzi`.

# Arguments:

  - `nxi`: Grid-points in `\\widehat{x}`-direction.

  - `nyi`: Grid-points in `\\widehat{y}`-direction.

  - `nzi`: Grid-points in `\\widehat{z}`-direction.
"""
struct TracerWKBImpact{A <: AbstractArray{<:AbstractFloat, 3}}
    uchi::A
    vchi::A
    wchi::A
    dchidt::A
end

function TracerWKBImpact(nxi::Integer, nyi::Integer, nzi::Integer)
    return TracerWKBImpact([zeros(nxi, nyi, nzi) for i in 1:4]...)
end
