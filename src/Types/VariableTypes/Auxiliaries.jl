"""
```julia
Auxiliaries{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractMatrix{<:AbstractFloat},
}
```

Container for the auxiliary array used in the reconstruction of prognostic variables and arrays used in the Thomas tridiagonal solver.

```julia
Auxiliaries(domain::Domain)::Auxiliaries
```

Construct an `Auxiliaries` instance with zero-initialized auxiliary arrays.

# Fields

  - `phi::A`: Auxiliary array used as input for [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref).

  - `ath::A`: Sub (lower) diagonal array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `bth::A`: Center diagonal array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `cth::A`: Super (upper) diagonal array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `fth::A`: Right-hand side array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `qth::A`: Work array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `pth::B`: Auxiliary array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `fth_bc::B`: Auxiliary right-hand side array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `qth_bc::B`: Auxiliary work array used as input for [`PinCFlow.Update.thomas_algorithm!`](@ref)

  - `stress_tensor_11::A`: Stress tensor element ``\\Pi^{1 1}``.

  - `stress_tensor_12::A`: Stress tensor element ``\\Pi^{1 2}``.
  
  - `stress_tensor_13::A`: Stress tensor element ``\\Pi^{1 3}``.
  
  - `stress_tensor_22::A`: Stress tensor element ``\\Pi^{2 2}``. 
  
  - `stress_tensor_23::A`: Stress tensor element ``\\Pi^{2 3}``. 
  
  - `stress_tensor_33::A`: Stress tensor element ``\\Pi^{3 3}``.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Auxiliaries{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractMatrix{<:AbstractFloat},
}
    phi::A
    ath::A
    bth::A
    cth::A
    fth::A
    qth::A
    pth::B
    fth_bc::B
    qth_bc::B
    stress_tensor_11::A
    stress_tensor_12::A 
    stress_tensor_13::A 
    stress_tensor_22::A 
    stress_tensor_23::A 
    stress_tensor_33::A 
end

function Auxiliaries(domain::Domain)::Auxiliaries
    (; nx, ny, nz, nxx, nyy, nzz) = domain

    return Auxiliaries(
        zeros(nxx, nyy, nzz),
        [zeros(nx, ny, nz) for i in 1:5]...,
        [zeros(nx, ny) for i in 1:3]...,
        [zeros(nxx, nyy, nzz) for i in 1:6]...,
    )
end
