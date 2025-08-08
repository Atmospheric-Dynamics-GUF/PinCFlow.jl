"""
```julia
Tensor{A <: AbstractArray{<:AbstractFloat, 3}}
```

Tensor elements of the linear operator, as computed by [`PinCFlow.PoissonSolver.compute_operator!`](@ref).

```julia
Tensor(domain::Domain)
```

Create a `Tensor` instance with zero-initialized arrays sized according to the dimensions of the MPI subdomain.

# Fields

- `ac_b::A`: Coefficient applied to ``s``.
- `acv_b::A`: Sum of vertical coefficients (upper and lower diagonals).
- `ach_b::A`: Sum of horizontal coefficients (left, right, forward and backward diagonals).
- `al_b::A`: Coefficient applied to ``s_{i - 1}``.
- `ar_b::A`: Coefficient applied to ``s_{i + 1}``.
- `ab_b::A`: Coefficient applied to ``s_{j - 1}``.
- `af_b::A`: Coefficient applied to ``s_{j + 1}``.
- `ad_b::A`: Coefficient applied to ``s_{k - 1}``.
- `au_b::A`: Coefficient applied to ``s_{k + 1}``.
- `ald_b::A`: Coefficient applied to ``s_{i - 1, k - 1}``.
- `alu_b::A`: Coefficient applied to ``s_{i - 1, k + 1}``.
- `ard_b::A`: Coefficient applied to ``s_{i + 1, k - 1}``.
- `aru_b::A`: Coefficient applied to ``s_{i + 1, k + 1}``.
- `abd_b::A`: Coefficient applied to ``s_{j - 1, k - 1}``.
- `abu_b::A`: Coefficient applied to ``s_{j - 1, k + 1}``.
- `afd_b::A`: Coefficient applied to ``s_{j + 1, k - 1}``.
- `afu_b::A`: Coefficient applied to ``s_{j + 1, k + 1}``.
- `add_b::A`: Coefficient applied to ``s_{k - 2}``.
- `auu_b::A`: Coefficient applied to ``s_{k + 2}``.
- `aldd_b::A`: Coefficient applied to ``s_{i - 1, k - 2}``.
- `aluu_b::A`: Coefficient applied to ``s_{i - 1, k + 2}``.
- `ardd_b::A`: Coefficient applied to ``s_{i + 1, k - 2}``.
- `aruu_b::A`: Coefficient applied to ``s_{i + 1, k + 2}``.
- `abdd_b::A`: Coefficient applied to ``s_{j - 1, k - 2}``.
- `abuu_b::A`: Coefficient applied to ``s_{j - 1, k + 2}``.
- `afdd_b::A`: Coefficient applied to ``s_{j + 1, k - 2}``.
- `afuu_b::A`: Coefficient applied to ``s_{j + 1, k + 2}``.

# Arguments

- `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Tensor{A <: AbstractArray{<:AbstractFloat, 3}}
    ac_b::A
    acv_b::A
    ach_b::A
    al_b::A
    ar_b::A
    ab_b::A
    af_b::A
    ad_b::A
    au_b::A
    ald_b::A
    alu_b::A
    ard_b::A
    aru_b::A
    abd_b::A
    abu_b::A
    afd_b::A
    afu_b::A
    add_b::A
    auu_b::A
    aldd_b::A
    aluu_b::A
    ardd_b::A
    aruu_b::A
    abdd_b::A
    abuu_b::A
    afdd_b::A
    afuu_b::A
end

function Tensor(domain::Domain)
    (; nx, ny, nz) = domain

    # Return a Tensor instance.
    return Tensor([zeros(nx, ny, nz) for i in 1:27]...)
end
