"""
```julia
BicGStab{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Workspace arrays for the BiCGStab iterative linear solver.

# Fields

  - `r_vm::A`: Vertically-averaged residual vector for convergence monitoring (nx × ny)
  - `p::B`: Search direction vector (nx × ny × nz)
  - `r0::B`: Initial residual vector (fixed throughout iteration) (nx × ny × nz)
  - `rold::B`: Previous residual vector (nx × ny × nz)
  - `r::B`: Current residual vector (nx × ny × nz)
  - `s::B`: Intermediate vector in BiCGStab algorithm (nx × ny × nz)
  - `t::B`: Matrix-vector product of s (nx × ny × nz)
  - `v::B`: Matrix-vector product of p (nx × ny × nz)
  - `matvec::B`: General matrix-vector product workspace (nx × ny × nz)
  - `v_pc::B`: Preconditioned vector workspace (nx × ny × nz)

# Usage

Provides temporary storage for [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref)
to avoid repeated memory allocations during iterative solving.
"""
struct BicGStab{
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

"""
```julia
BicGStab(namelists::Namelists, domain::Domain)
```

Initialize BiCGStab workspace arrays sized according to local domain.

# Arguments

  - `namelists::Namelists`: Configuration parameters (unused)
  - `domain::Domain`: Local domain dimensions

# Returns

  - `BicGStab`: Workspace container with zero-initialized arrays
"""
function BicGStab(namelists::Namelists, domain::Domain)
    (; nx, ny, nz) = domain
    return BicGStab(zeros(nx, ny), [zeros(nx, ny, nz) for i in 1:9]...)
end
