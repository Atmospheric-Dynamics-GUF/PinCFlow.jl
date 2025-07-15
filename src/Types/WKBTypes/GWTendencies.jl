"""
```julia
GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
```

Gravity wave drag and heating tendency fields.

# Fields

  - `dudt`: Zonal wind tendency (gravity wave drag)
  - `dvdt`: Meridional wind tendency (gravity wave drag)
  - `dthetadt`: Potential temperature tendency (gravity wave heating)
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

Construct a `GWTendencies` instance.

# Arguments

  - `nxx`: Number of grid points in x-direction
  - `nyy`: Number of grid points in y-direction
  - `nzz`: Number of grid points in z-direction

# Returns

  - `GWTendencies`: Initialized container with zero arrays for gravity wave tendencies
"""
function GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)
    return GWTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
