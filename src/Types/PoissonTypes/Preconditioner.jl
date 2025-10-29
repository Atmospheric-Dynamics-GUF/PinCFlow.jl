"""
```julia
Preconditioner{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractMatrix{<:AbstractFloat},
    E <: AbstractMatrix{<:AbstractFloat},
}
```

Workspace arrays for applying the preconditioner.

```julia
Preconditioner(domain::Domain)::Preconditioner
```

Create a `Preconditioner` instance with zero-initialized arrays sized according to the dimensions of the MPI subdomain.

# Fields

  - `s_pc::A`: Solution computed by the preconditioner.

  - `q_pc::B`: Auxiliary array used for the upward sweep.

  - `p_pc::C`: Auxiliary array used for the upward sweep and downward pass.

  - `s_pc_bc::D`: MPI communication buffer for `s_pc`.

  - `q_pc_bc::E`: MPI communication buffer for `q_pc`.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Preconditioner{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractMatrix{<:AbstractFloat},
    E <: AbstractMatrix{<:AbstractFloat},
}
    s_pc::A
    q_pc::B
    p_pc::C
    s_pc_bc::D
    q_pc_bc::E
end

function Preconditioner(domain::Domain)::Preconditioner
    (; nx, ny, nz) = domain

    return Preconditioner(
        [zeros(nx, ny, nz) for i in 1:2]...,
        [zeros(nx, ny) for i in 1:3]...,
    )
end
