"""
```julia
WKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Integrals of ray-volume properties.

```julia
WKBIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)::WKBIntegrals
```

Construct a `WKBIntegrals` instance, with arrays sized according to the given dimensions.

# Fields

  - `uu::A`: Zonal zonal-momentum flux.

  - `uv::A`: Meridional zonal-momentum flux.

  - `uw::A`: Vertical zonal-momentum flux.

  - `vv::A`: Meridional meridional-momentum flux.

  - `vw::A`: Vertical meridional-momentum flux.

  - `utheta::A`: Zonal mass-weighted potential-temperature flux.

  - `vtheta::A`: Meridional mass-weighted potential-temperature flux.

  - `e::A`: Wave-energy density.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct WKBIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
    uu::A
    uv::A
    uw::A
    vv::A
    vw::A
    utheta::A
    vtheta::A
    e::A
    sterm::A
    bterm::A
    q00::A
    q10::A
    q20::A
end

function WKBIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)::WKBIntegrals
    return WKBIntegrals([zeros(nxx, nyy, nzz) for i in 1:13]...)
end
