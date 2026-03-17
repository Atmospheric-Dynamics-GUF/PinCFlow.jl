"""
```julia
Auxiliaries{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractMatrix{<:AbstractFloat},
}
```

Container for the auxiliary array used in the reconstruction of prognostic variables and arrays used in the Thomas tridiagonal solver.

```julia
Auxiliaries(domain::Domain)::Auxiliaries
```

Construct an `Auxiliaries` instance with zero-initialized auxiliary arrays.

# Fields

  - `phi::A`: Auxiliary array used as input for [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref).

  - `ath::B`: Sub (lower) diagonal array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `bth::B`: Center diagonal array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `cth::B`: Super (upper) diagonal array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `fth::B`: Right-hand side array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `qth::B`: Work array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `pth::C`: Auxiliary array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `fth_bc::C`: Auxiliary right-hand side array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `qth_bc::C`: Auxiliary work array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

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
