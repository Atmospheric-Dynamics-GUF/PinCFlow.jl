"""
```julia
Preconditioner{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractMatrix{<:AbstractFloat},
}
```

Workspace arrays for line relaxation preconditioner using alternating direction implicit (ADI) method.

# Fields

  - `s_pc::A`: Preconditioned solution workspace (nx × ny × nz)
  - `q_pc::A`: Upper diagonal coefficients for tridiagonal solve (nx × ny × nz)
  - `p_pc::B`: Diagonal pivot workspace for elimination (nx × ny)
  - `s_pc_bc::B`: Boundary communication buffer for solution (nx × ny)
  - `q_pc_bc::B`: Boundary communication buffer for coefficients (nx × ny)

# Usage

Provides temporary storage for [`PinCFlow.PoissonSolver.apply_preconditioner!`](@ref)
to perform vertical line relaxation with MPI communication between domains.
"""
struct Preconditioner{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractMatrix{<:AbstractFloat},
}
    s_pc::A
    q_pc::A
    p_pc::B
    s_pc_bc::B
    q_pc_bc::B
end

"""
```julia
Preconditioner(domain::Domain)
```

Initialize preconditioner workspace arrays sized according to local domain.

# Arguments

  - `domain::Domain`: Local domain dimensions

# Returns

  - `Preconditioner`: Container with zero-initialized workspace arrays
"""
function Preconditioner(domain::Domain)

    # Get all necessary fields.
    (; nx, ny, nz) = domain

    # Return a Preconditioner instance.
    return Preconditioner(
        [zeros(nx, ny, nz) for i in 1:2]...,
        [zeros(nx, ny) for i in 1:3]...,
    )
end
