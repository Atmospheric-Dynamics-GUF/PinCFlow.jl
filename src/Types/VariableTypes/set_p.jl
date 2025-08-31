"""
```julia 
set_p(
    model::AbstractModel,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pstrattfc::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}

Construct `p` with zero-size array for non-compressible dynamic equations.

```julia 
set_p(
    model::Compressible,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pstrattfc::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}

Construct `p` with as copy of `pstrattfc` in compressible mode.

```
"""

function set_p end

function set_p(
    model::AbstractModel,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pstrattfc::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return zeros(0, 0, 0)
end

function set_p(
    model::Compressible,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pstrattfc::AbstractArray{<:AbstractFloat, 3},
)::AbstractArray{<:AbstractFloat, 3}
    return copy(pstrattfc)
end