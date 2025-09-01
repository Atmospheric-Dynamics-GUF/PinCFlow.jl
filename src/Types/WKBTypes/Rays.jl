"""
```julia
Rays{A <: AbstractArray{<:AbstractFloat, 4}}
```

Container for prognostic ray-volume properties.

```julia
Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
```

Construct a `Rays` instance, with arrays sized according to the given dimensions.

# Fields

  - `x::A`: Position in ``x``.
  - `y::A`: Position in ``y``.
  - `z::A`: Position in ``z``.
  - `k::A`: Position in ``k``.
  - `l::A`: Position in ``l``.
  - `m::A`: Position in ``m``.
  - `dxray::A`: Extent in ``x``.
  - `dyray::A`: Extent in ``y``.
  - `dzray::A`: Extent in ``z``.
  - `dkray::A`: Extent in ``k``.
  - `dlray::A`: Extent in ``l``.
  - `dmray::A`: Extent in ``m``.
  - `dens::A`: Phase-space wave-action density.
  - `dphi::A`: phase.

# Arguments

  - `nray_wrk`: Size of the spectral dimension of ray-volume arrays.
  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.
  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct Rays{A <: AbstractArray{<:AbstractFloat, 4}}
    x::A
    y::A
    z::A
    k::A
    l::A
    m::A
    dxray::A
    dyray::A
    dzray::A
    dkray::A
    dlray::A
    dmray::A
    dens::A
    dphi::A
end

function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer, namelists::Namelists)
    (; icesetup) = namelists.ice
    return Rays(nray_wrk, nxx, nyy, nzz, icesetup)
end

function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer, icesetup::NoIce)
    return Rays([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:14]...)
end

function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer, icesetup::AbstractIce)
    return Rays([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:14]...)
end