"""
```julia
WKBAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Auxiliary for calculating the gravity-wave induced shear in the turbulence scheme.

```julia
WKBAuxiliaries(nxx::Integer, nyy::Integer, nzz::Integer)::WKBIntegrals
```

Construct a `WKBAuxiliaries` instance, with an array sized according to the given dimensions.

# Fields

  - `shear::A`: Gravity-wave induced shear field.

# Arguments

  - `nxx`: Number of subdomain grid points in ``\\hat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\hat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\hat{z}``-direction.
"""
struct WKBAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    shear::A
end

function WKBAuxiliaries(nxx::Integer, nyy::Integer, nzz::Integer)::WKBAuxiliaries
    return WKBAuxiliaries(zeros(nxx, nyy, nzz))
end
