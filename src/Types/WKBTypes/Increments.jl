"""
```julia
Increments{A <: AbstractArray{<:AbstractFloat, 4}}
```

Ray-volume-propagation increments.

```julia
Increments(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
```

Construct an `Increments` instance, with arrays sized according to the given dimensions.

# Fields

  - `dxray::A`: Increments for the position in ``x``.
  - `dyray::A`: Increments for the position in ``y``.
  - `dzray::A`: Increments for the position in ``z``.
  - `dkray::A`: Increments for the position in ``k``.
  - `dlray::A`: Increments for the position in ``l``.
  - `dmray::A`: Increments for the position in ``m``.
  - `ddxray::A`: Increments for the extent in ``x``.
  - `ddyray::A`: Increments for the extent in ``y``.
  - `ddzray::A`: Increments for the extent in ``z``.

# Arguments

  - `nray_wrk`: Size of the spectral dimension of ray-volume arrays.
  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.
  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
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

function Increments(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
    return Increments([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:9]...)
end
