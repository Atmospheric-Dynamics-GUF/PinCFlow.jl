"""
```julia
Increments{A <: AbstractArray{<:AbstractFloat, 4}}
```

Ray propagation increments for position and wavenumber evolution.

# Fields

  - `dxray`, `dyray`, `dzray`: Position increments
  - `dkray`, `dlray`, `dmray`: Wavenumber increments
  - `ddxray`, `ddyray`, `ddzray`: Second-order position increments
"""
struct Increments{A <: AbstractArray{<:AbstractFloat, 4}}
    dxray::A
    dyray::A
    dzray::A
    dkray::A
    dlray::A
    dmray::A
    ddxray::A
    ddyray::A
    ddzray::A
end

"""
```julia
Increments(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
```

Construct an `Increments` instance.
"""
function Increments(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
    return Increments([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:9]...)
end
