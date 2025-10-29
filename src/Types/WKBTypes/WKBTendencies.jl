"""
```julia
WKBTendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
```

Gravity-wave drag and heating fields.

```julia
WKBTendencies(nxx::Integer, nyy::Integer, nzz::Integer)::WKBTendencies
```

Construct a `WKBTendencies` instance, with arrays sized according to the given dimensions.

# Fields

  - `dudt::A`: Gravity-wave drag on the zonal momentum.

  - `dvdt::B`: Gravity-wave drag on the meridional momentum.

  - `dthetadt::C`: Gravity-wave heating term in the mass-weighted potential-temperature equation.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\widehat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\widehat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\widehat{z}``-direction.
"""
struct WKBTendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
}
    dudt::A
    dvdt::B
    dthetadt::C
end

function WKBTendencies(nxx::Integer, nyy::Integer, nzz::Integer)::WKBTendencies
    return WKBTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
