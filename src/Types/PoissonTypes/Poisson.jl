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

Main container for Poisson-solver workspace and solution arrays.

# Fields

  - `rhs::A`: Right-hand side.
  - `solution::A`: Solution of the Poisson problem.
  - `tensor::B`: Tensor elements of the linear operator.
  - `operator::C`: Workspace arrays for applying the linear operator.
  - `preconditioner::D`: Workspace arrays for applying the preconditioner.
  - `bicgstab::E`: Workspace arrays used by the BicGStab algorithm.
  - `correction::F`: Correction terms used to update the horizontal wind in the corrector step.
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
Poisson(domain::Domain)
```

Initialize the complete Poisson-solver workspace, sized according to the dimensions of the MPI subdomain.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Returns

  - `::Poisson`: `Poisson` instance with zero-initialized arrays.

# See also

  - [`PinCFlow.Types.PoissonTypes.Tensor`](@ref)
  - [`PinCFlow.Types.PoissonTypes.Operator`](@ref)
  - [`PinCFlow.Types.PoissonTypes.Preconditioner`](@ref)
  - [`PinCFlow.Types.PoissonTypes.BicGStab`](@ref)
  - [`PinCFlow.Types.PoissonTypes.Correction`](@ref)
"""
function Poisson(domain::Domain)

    # Get all necessary fields.
    (; nx, ny, nz) = domain

    # Initialize everything.
    (rhs, solution) = (zeros(nx, ny, nz) for i in 1:2)
    tensor = Tensor(domain)
    operator = Operator(domain)
    preconditioner = Preconditioner(domain)
    bicgstab = BicGStab(domain)
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
