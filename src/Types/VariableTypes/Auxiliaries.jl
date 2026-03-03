"""
```julia
Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Auxiliary array used in the reconstruction of prognostic variables.

```julia
Auxiliaries(domain::Domain)::Auxiliaries
```

Construct an `Auxiliaries` instance with a zero-initialized auxiliary array sized according to the MPI subdomain dimensions.

# Fields

  - `phi::A`: Auxiliary array used as input for [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref).

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Auxiliaries{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractMatrix{<:AbstractFloat},
}
    phi::A
    ath::B
    bth::B
    cth::B
    fth::B
    qth::B
    pth::C
    fth_bc::C
    qth_bc::C
end

function Auxiliaries(domain::Domain)::Auxiliaries
    (; nx, ny, nz, nxx, nyy, nzz) = domain

    return Auxiliaries(
        zeros(nxx, nyy, nzz),
        [zeros(nx, ny, nz) for i in 1:5]...,
        [zeros(nx, ny) for i in 1:3]...,
    )
end
