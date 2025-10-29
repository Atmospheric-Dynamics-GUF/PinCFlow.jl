"""
```julia
BicGStab{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 3},
    H <: AbstractArray{<:AbstractFloat, 3},
    I <: AbstractArray{<:AbstractFloat, 3},
    J <: AbstractArray{<:AbstractFloat, 3},
}
```

Workspace arrays used by [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref).

```julia
BicGStab(domain::Domain)::BicGStab
```

Create a `BicGStab` instance with zero-initialized workspace arrays sized according to dimensions of the MPI subdomain.

# Fields

  - `r_vm::A`: Vertically-averaged residual.

  - `p::B`: Search direction.

  - `r0::C`: Initial residual.

  - `rold::D`: Previous residual.

  - `r::E`: Current residual.

  - `s::F`: Intermediate solution.

  - `t::G`: Result of applying the linear operator to `s`.

  - `v::H`: Result of applying the linear operator to `p`.

  - `matvec::I`: Intermediate result of applying the linear operator.

  - `v_pc::J`: Output of the preconditioner.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct BicGStab{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractArray{<:AbstractFloat, 3},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 3},
    G <: AbstractArray{<:AbstractFloat, 3},
    H <: AbstractArray{<:AbstractFloat, 3},
    I <: AbstractArray{<:AbstractFloat, 3},
    J <: AbstractArray{<:AbstractFloat, 3},
}
    r_vm::A
    p::B
    r0::C
    rold::D
    r::E
    s::F
    t::G
    v::H
    matvec::I
    v_pc::J
end

function BicGStab(domain::Domain)::BicGStab
    (; nx, ny, nz) = domain

    return BicGStab(zeros(nx, ny), [zeros(nx, ny, nz) for i in 1:9]...)
end
