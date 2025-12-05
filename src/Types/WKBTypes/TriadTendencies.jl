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
struct TriadTendencies{A <: AbstractArray{<: AbstractFloat, 5}}
    wavespectrum::A  
end

function TriadTendencies(nxx::Integer, nyy::Integer, nzz::Integer, kp_size::Integer, m_size::Integer)::TriadTendencies
    return TriadTendencies(zeros(nxx, nyy, nzz, kp_size, m_size))
end
