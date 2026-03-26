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

  - `nxx`: Number of subdomain grid points in ``\\hat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\hat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\hat{z}``-direction.
"""
struct WKBIntegrals{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:ComplexF64, 3},
}
    uu::A
    uv::A
    uw::A
    vv::A
    vw::A
    utheta::A
    vtheta::A
    e::A
    uhat::B
    vhat::B
    what::B
    bhat::B
    pihat::B
    chihat::B
end

function WKBIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)::WKBIntegrals
    return WKBIntegrals(
        [zeros(nxx, nyy, nzz) for i in 1:8]...,
        [zeros(ComplexF64, nxx, nyy, nzz) for i in 1:6]...,
    )
end
