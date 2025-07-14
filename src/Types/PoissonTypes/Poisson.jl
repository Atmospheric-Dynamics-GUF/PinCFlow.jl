"""
```julia
Poisson{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: Tensor,
    C <: Operator,
    D <: Preconditioner,
    E <: BicGStab,
    F <: Correction,
}
```

Main container for Poisson solver workspace and solution arrays.

# Fields

  - `rhs::A`: Right-hand side vector
  - `solution::A`: Solution vector
  - `tensor::B`: Matrix coefficient tensor for 27-point stencil
  - `operator::C`: Operator workspace for matrix-vector products
  - `preconditioner::D`: Preconditioner workspace arrays
  - `bicgstab::E`: BiCGStab iterative solver workspace
  - `correction::F`: Pressure correction terms for velocity field

# Usage

Central data structure used by [`PinCFlow.PoissonSolver.solve_poisson!`](@ref)
to solve the pressure Poisson equation and apply corrections to maintain
mass conservation.
"""
struct Poisson{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: Tensor,
    C <: Operator,
    D <: Preconditioner,
    E <: BicGStab,
    F <: Correction,
}
    rhs::A
    solution::A
    tensor::B
    operator::C
    preconditioner::D
    bicgstab::E
    correction::F
end

"""
```julia
Poisson(namelists::Namelists, domain::Domain)
```

Initialize complete Poisson solver workspace sized according to local domain.

# Arguments

  - `namelists::Namelists`: Configuration parameters
  - `domain::Domain`: Local domain dimensions

# Returns

  - `Poisson`: Container with all solver components initialized
"""
function Poisson(namelists::Namelists, domain::Domain)

    # Get all necessary fields.
    (; nx, ny, nz) = domain

    # Initialize everything.
    (rhs, solution) = (zeros(nx, ny, nz) for i in 1:2)
    tensor = Tensor(domain)
    operator = Operator(domain)
    preconditioner = Preconditioner(domain)
    bicgstab = BicGStab(namelists, domain)
    correction = Correction(domain)

    # Return a Poisson instance.
    return Poisson(
        rhs,
        solution,
        tensor,
        operator,
        preconditioner,
        bicgstab,
        correction,
    )
end
