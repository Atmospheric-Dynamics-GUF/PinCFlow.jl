"""
```julia
BiCGSTAB{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Workspace arrays used by [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref).

```julia
BiCGSTAB(domain::Domain)::BiCGSTAB
```

Create a `BiCGSTAB` instance with zero-initialized workspace arrays sized according to dimensions of the MPI subdomain.

# Fields

  - `r_vm::A`: Vertically-averaged residual.

  - `p::B`: Search direction.

  - `r0::B`: Initial residual.

  - `rold::B`: Previous residual.

  - `r::B`: Current residual.

  - `s::B`: Intermediate solution.

  - `t::B`: Result of applying the linear operator to `s`.

  - `v::B`: Result of applying the linear operator to `p`.

  - `matvec::B`: Intermediate result of applying the linear operator.

  - `v_pc::B`: Output of the preconditioner.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct BiCGSTAB{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    r_vm::A
    p::B
    r0::B
    rold::B
    r::B
    s::B
    t::B
    v::B
    matvec::B
    v_pc::B
end

function BiCGSTAB(domain::Domain)::BiCGSTAB
    (; nx, ny, nz) = domain

    return BiCGSTAB(zeros(nx, ny), [zeros(nx, ny, nz) for i in 1:9]...)
end
