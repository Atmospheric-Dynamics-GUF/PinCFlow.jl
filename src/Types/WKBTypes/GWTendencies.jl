"""
```julia
GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Gravity-wave drag and heating fields.

# Fields

  - `dudt::A`: Gravity-wave drag on the zonal momentum.
  - `dvdt::A`: Gravity-wave drag on the meridional momentum.
  - `dthetadt::A`: Gravity-wave heating term in the mass-weighted potential-temperature equation.
"""
struct GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dudt::A
    dvdt::A
    dthetadt::A
end

"""
```julia
GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)
```

Construct a `GWTendencies` instance, with arrays sized according to the given dimensions.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.
  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.
  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.

# Returns

  - `::GWTendencies`: `GWTendencies` instance with zero-initialized arrays.
"""
function GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)
    return GWTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
