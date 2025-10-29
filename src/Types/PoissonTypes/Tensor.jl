"""
```julia
Tensor{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 3},
    H <: AbstractArray{<:AbstractFloat, 3},
    I <: AbstractArray{<:AbstractFloat, 3},
    J <: AbstractArray{<:AbstractFloat, 3},
    K <: AbstractArray{<:AbstractFloat, 3},
    L <: AbstractArray{<:AbstractFloat, 3},
    M <: AbstractArray{<:AbstractFloat, 3},
    N <: AbstractArray{<:AbstractFloat, 3},
    O <: AbstractArray{<:AbstractFloat, 3},
    P <: AbstractArray{<:AbstractFloat, 3},
    Q <: AbstractArray{<:AbstractFloat, 3},
    R <: AbstractArray{<:AbstractFloat, 3},
    S <: AbstractArray{<:AbstractFloat, 3},
    T <: AbstractArray{<:AbstractFloat, 3},
    U <: AbstractArray{<:AbstractFloat, 3},
    V <: AbstractArray{<:AbstractFloat, 3},
    W <: AbstractArray{<:AbstractFloat, 3},
    X <: AbstractArray{<:AbstractFloat, 3},
    Y <: AbstractArray{<:AbstractFloat, 3},
}
```

Tensor elements of the linear operator, as computed by [`PinCFlow.PoissonSolver.compute_operator!`](@ref).

```julia
Tensor(domain::Domain)::Tensor
```

Create a `Tensor` instance with zero-initialized arrays sized according to the dimensions of the MPI subdomain.

# Fields

  - `ac_b::A`: Coefficient applied to ``s``.

  - `al_b::B`: Coefficient applied to ``s_{i - 1}``.

  - `ar_b::C`: Coefficient applied to ``s_{i + 1}``.

  - `ab_b::D`: Coefficient applied to ``s_{j - 1}``.

  - `af_b::E`: Coefficient applied to ``s_{j + 1}``.

  - `ad_b::F`: Coefficient applied to ``s_{k - 1}``.

  - `au_b::G`: Coefficient applied to ``s_{k + 1}``.

  - `ald_b::H`: Coefficient applied to ``s_{i - 1, k - 1}``.

  - `alu_b::I`: Coefficient applied to ``s_{i - 1, k + 1}``.

  - `ard_b::J`: Coefficient applied to ``s_{i + 1, k - 1}``.

  - `aru_b::K`: Coefficient applied to ``s_{i + 1, k + 1}``.

  - `abd_b::L`: Coefficient applied to ``s_{j - 1, k - 1}``.

  - `abu_b::M`: Coefficient applied to ``s_{j - 1, k + 1}``.

  - `afd_b::N`: Coefficient applied to ``s_{j + 1, k - 1}``.

  - `afu_b::O`: Coefficient applied to ``s_{j + 1, k + 1}``.

  - `add_b::P`: Coefficient applied to ``s_{k - 2}``.

  - `auu_b::Q`: Coefficient applied to ``s_{k + 2}``.

  - `aldd_b::R`: Coefficient applied to ``s_{i - 1, k - 2}``.

  - `aluu_b::S`: Coefficient applied to ``s_{i - 1, k + 2}``.

  - `ardd_b::T`: Coefficient applied to ``s_{i + 1, k - 2}``.

  - `aruu_b::U`: Coefficient applied to ``s_{i + 1, k + 2}``.

  - `abdd_b::V`: Coefficient applied to ``s_{j - 1, k - 2}``.

  - `abuu_b::W`: Coefficient applied to ``s_{j - 1, k + 2}``.

  - `afdd_b::X`: Coefficient applied to ``s_{j + 1, k - 2}``.

  - `afuu_b::Y`: Coefficient applied to ``s_{j + 1, k + 2}``.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Tensor{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 3},
    H <: AbstractArray{<:AbstractFloat, 3},
    I <: AbstractArray{<:AbstractFloat, 3},
    J <: AbstractArray{<:AbstractFloat, 3},
    K <: AbstractArray{<:AbstractFloat, 3},
    L <: AbstractArray{<:AbstractFloat, 3},
    M <: AbstractArray{<:AbstractFloat, 3},
    N <: AbstractArray{<:AbstractFloat, 3},
    O <: AbstractArray{<:AbstractFloat, 3},
    P <: AbstractArray{<:AbstractFloat, 3},
    Q <: AbstractArray{<:AbstractFloat, 3},
    R <: AbstractArray{<:AbstractFloat, 3},
    S <: AbstractArray{<:AbstractFloat, 3},
    T <: AbstractArray{<:AbstractFloat, 3},
    U <: AbstractArray{<:AbstractFloat, 3},
    V <: AbstractArray{<:AbstractFloat, 3},
    W <: AbstractArray{<:AbstractFloat, 3},
    X <: AbstractArray{<:AbstractFloat, 3},
    Y <: AbstractArray{<:AbstractFloat, 3},
}
    ac_b::A
    al_b::B
    ar_b::C
    ab_b::D
    af_b::E
    ad_b::F
    au_b::G
    ald_b::H
    alu_b::I
    ard_b::J
    aru_b::K
    abd_b::L
    abu_b::M
    afd_b::N
    afu_b::O
    add_b::P
    auu_b::Q
    aldd_b::R
    aluu_b::S
    ardd_b::T
    aruu_b::U
    abdd_b::V
    abuu_b::W
    afdd_b::X
    afuu_b::Y
end

function Tensor(domain::Domain)::Tensor
    (; nx, ny, nz) = domain

    return Tensor([zeros(nx, ny, nz) for i in 1:25]...)
end
