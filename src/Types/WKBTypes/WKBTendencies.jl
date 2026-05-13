"""
```julia
WKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Gravity-wave drag and heating fields.

```julia
WKBTendencies(nxx::Integer, nyy::Integer, nzz::Integer)::WKBTendencies
```

Construct a `WKBTendencies` instance, with arrays sized according to the given dimensions.

# Fields

  - `dudt::A`: Gravity-wave drag on the zonal momentum.

  - `dvdt::A`: Gravity-wave drag on the meridional momentum.

  - `dthetadt::A`: Gravity-wave heating term in the mass-weighted potential-temperature equation.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct WKBTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dudt::A
    dvdt::A
    dthetadt::A
    dtkedt::A
end

function WKBTendencies(nxx::Integer, nyy::Integer, nzz::Integer)::WKBTendencies
    return WKBTendencies([zeros(nxx, nyy, nzz) for i in 1:4]...)
end
