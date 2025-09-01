"""
```julia
WKBIncrements{A <: AbstractArray{<:AbstractFloat, 4}}
```

Ray-volume-propagation increments.

```julia
WKBIncrements(
    nray_wrk::Integer,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::WKBIncrements
```

Construct an `WKBIncrements` instance, with arrays sized according to the given dimensions.

# Fields

  - `dxray::A`: WKBIncrements for the position in ``x``.

  - `dyray::A`: WKBIncrements for the position in ``y``.

  - `dzray::A`: WKBIncrements for the position in ``z``.

  - `dkray::A`: WKBIncrements for the position in ``k``.

  - `dlray::A`: WKBIncrements for the position in ``l``.

  - `dmray::A`: WKBIncrements for the position in ``m``.

  - `ddxray::A`: WKBIncrements for the extent in ``x``.

  - `ddyray::A`: WKBIncrements for the extent in ``y``.

  - `ddzray::A`: WKBIncrements for the extent in ``z``.

  - `dpray::A`: Increments for the phase ``dphi``.


# Arguments

  - `nray_wrk`: Size of the spectral dimension of ray-volume arrays.

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct WKBIncrements{A <: AbstractArray{<:AbstractFloat, 4}}
    dxray::A
    dyray::A
    dzray::A
    dkray::A
    dlray::A
    dmray::A
    ddxray::A
    ddyray::A
    ddzray::A
    dpray::A
end

function WKBIncrements(
    nray_wrk::Integer,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::WKBIncrements
    return WKBIncrements([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:10]...)
end
