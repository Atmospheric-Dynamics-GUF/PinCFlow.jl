"""
```julia
GWIntegrals{A <: AbstractArray{<:AbstractFloat, 3}}
```

Integrals of ray-volume properties.

# Fields

  - `uu::A`: Zonal zonal-momentum flux.
  - `uv::A`: Meridional zonal-momentum flux.
  - `uw::A`: Vertical zonal-momentum flux.
  - `vv::A`: Meridional meridional-momentum flux.
  - `vw::A`: Vertical meridional-momentum flux.
  - `etx::A`: Elastic term in the zonal-momentum equation.
  - `ety::A`: Elastic term in the meridional momentum equation.
  - `utheta::A`: Zonal mass-weighted potential-temperature flux.
  - `vtheta::A`: Meridional mass-weighted potential-temperature flux.
  - `e::A`: Wave-energy density.
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

Construct a `GWIntegrals` instance, with arrays sized according to the given dimensions.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.
  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.

# Returns

  - `::GWIntegrals`: `GWIntegrals` instance with zero-initialized arrays.
"""
function GWIntegrals(nxx::Integer, nyy::Integer, nzz::Integer)
    return GWIntegrals([zeros(nxx, nyy, nzz) for i in 1:10]...)
end
