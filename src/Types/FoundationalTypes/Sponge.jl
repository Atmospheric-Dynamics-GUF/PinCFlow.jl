"""
```julia
Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
```

Composite type for Rayleigh-damping coefficients and an auxiliary array for the computation of horizontal means.

```julia
Sponge(namelists::Namelists, domain::Domain)::Sponge
```

Construct a `Sponge` instance with zero-initialized arrays.

# Fields

  - `alphar::A`: Coefficient of the LHS sponge (used in all prognostic equations).

  - `betar::A`: Coefficient of the RHS sponge (used in the momentum equation).

  - `horizontal_mean::C`: Auxiliary array for the computation of horizontal means.

# Arguments

  - `namelists`: Namelists with all model parameters.

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

function Sponge(namelists::Namelists, domain::Domain)::Sponge
    (; float_type) = namelists.discretization
    (; nxx, nyy, nzz, nz) = domain

    return Sponge(
        [zeros(float_type, nxx, nyy, nzz) for i in 1:2]...,
        zeros(float_type, nz),
    )
end
