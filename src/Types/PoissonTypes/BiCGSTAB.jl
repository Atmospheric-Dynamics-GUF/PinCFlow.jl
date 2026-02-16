"""
```julia
BiCGSTAB{
    A <: AbstractArray{<:AbstractFloat, 3},
}
```

Workspace arrays used by [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref).

```julia
BiCGSTAB(domain::Domain)::BiCGSTAB
```

Create a `BiCGSTAB` instance with zero-initialized workspace arrays sized according to dimensions of the MPI subdomain.

# Fields

  - `p::A`: Search direction.

  - `r0::A`: Initial residual.

  - `rold::A`: Previous residual.

  - `r::A`: Current residual.

  - `s::A`: Intermediate solution.

  - `t::A`: Result of applying the linear operator to `s`.

  - `v::A`: Result of applying the linear operator to `p`.

  - `matvec::A`: Intermediate result of applying the linear operator.

  - `v_pc::A`: Output of the preconditioner.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct BiCGSTAB{
    A <: AbstractArray{<:AbstractFloat, 3},
}
    p::A
    r0::A
    rold::A
    r::A
    s::A
    t::A
    v::A
    matvec::A
    v_pc::A
end

function BiCGSTAB(domain::Domain)::BiCGSTAB
    (; nx, ny, nz) = domain

    return BiCGSTAB([zeros(nx, ny, nz) for i in 1:9]...)
end
