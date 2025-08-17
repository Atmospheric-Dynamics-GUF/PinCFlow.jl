"""
```julia
GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Gravity-wave drag and heating fields.

```julia
GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)::GWTendencies
```

Construct a `GWTendencies` instance, with arrays sized according to the given dimensions.

# Fields

  - `dudt::A`: Gravity-wave drag on the zonal momentum.

  - `dvdt::A`: Gravity-wave drag on the meridional momentum.

  - `dthetadt::A`: Gravity-wave heating term in the mass-weighted potential-temperature equation.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dudt::A
    dvdt::A
    dthetadt::A
end

function GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)::GWTendencies
    return GWTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
