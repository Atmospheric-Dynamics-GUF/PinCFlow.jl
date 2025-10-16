"""
```julia
Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
```

Composite type for Rayleigh-damping coefficients and an auxiliary array for the computation of horizontal means.

```julia
Sponge(domain::Domain)::Sponge
```

Construct a `Sponge` instance with zero-initialized arrays.

# Fields

  - `alphar::A`: Coefficient of the LHS sponge (used in all prognostic equations).

  - `betar::A`: Coefficient of the RHS sponge (used in the momentum equation).

  - `horizontal_mean::C`: Auxiliary array for the computation of horizontal means.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
    alphar::A
    betar::A
    horizontal_mean::B
end

function Sponge(domain::Domain)::Sponge
    (; nxx, nyy, nzz, nz) = domain

    return Sponge([zeros(nxx, nyy, nzz) for i in 1:2]..., zeros(nz))
end
