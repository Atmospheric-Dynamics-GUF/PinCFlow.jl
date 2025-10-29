"""
```julia
WKBIntegrals{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 3},
    H <: AbstractArray{<:AbstractFloat, 3},
}
```

Integrals of ray-volume properties.

```julia
WKBIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)::WKBIntegrals
```

Construct a `WKBIntegrals` instance, with arrays sized according to the given dimensions.

# Fields

  - `uu::A`: Zonal zonal-momentum flux.

  - `uv::B`: Meridional zonal-momentum flux.

  - `uw::C`: Vertical zonal-momentum flux.

  - `vv::D`: Meridional meridional-momentum flux.

  - `vw::E`: Vertical meridional-momentum flux.

  - `utheta::F`: Zonal mass-weighted potential-temperature flux.

  - `vtheta::G`: Meridional mass-weighted potential-temperature flux.

  - `e::H`: Wave-energy density.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct WKBIntegrals{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 3},
    H <: AbstractArray{<:AbstractFloat, 3},
}
    uu::A
    uv::B
    uw::C
    vv::D
    vw::E
    utheta::F
    vtheta::G
    e::H
end

function WKBIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)::WKBIntegrals
    return WKBIntegrals([zeros(nxx, nyy, nzz) for i in 1:8]...)
end
