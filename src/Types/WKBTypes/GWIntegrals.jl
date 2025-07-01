"""
```julia
GWIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Gravity wave momentum and energy integral quantities.

# Fields

  - `uu`, `uv`, `uw`, `vv`, `vw`: Momentum flux components
  - `etx`, `ety`: Horizontal energy flux components
  - `utheta`, `vtheta`: Temperature flux components
  - `e`: Wave energy density
"""
struct GWIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    uu::A
    uv::A
    uw::A
    vv::A
    vw::A
    etx::A
    ety::A
    utheta::A
    vtheta::A
    e::A
end

"""
```julia
GWIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)
```

Construct a `GWIntegrals` instance.
"""
function GWIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)
    return GWIntegrals([zeros(nxx, nyy, nzz) for i in 1:10]...)
end
