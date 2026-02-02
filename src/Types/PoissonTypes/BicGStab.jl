"""
```julia
BicGStab{
    A <: AbstractArray{<:AbstractFloat, 3},
}
```

Workspace arrays used by [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref).

```julia
BicGStab(domain::Domain)::BicGStab
```

Create a `BicGStab` instance with zero-initialized workspace arrays sized according to dimensions of the MPI subdomain.

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
struct BicGStab{
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

function BicGStab(domain::Domain)::BicGStab
    (; nx, ny, nz) = domain

    return BicGStab([zeros(nx, ny, nz) for i in 1:9]...)
end
