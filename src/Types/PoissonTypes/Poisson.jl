"""
```julia
Poisson{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: Tensor,
    C <: Operator,
    D <: BiCGSTAB,
    E <: Correction,
}
```

Main container for Poisson-solver workspace and solution arrays.

```julia
Poisson(domain::Domain)::Poisson
```

Create a `Poisson` instance with an initialized Poisson-solver workspace, sized according to the dimensions of the MPI subdomain.

# Fields

  - `rhs::A`: Right-hand side.

  - `solution::A`: Solution of the Poisson problem.

  - `tensor::B`: Tensor elements of the linear operator.

  - `operator::C`: Workspace arrays for applying the linear operator.

  - `bicgstab::D`: Workspace arrays used by the BiCGSTAB algorithm.

  - `correction::E`: Correction terms used to update the horizontal wind in the corrector step.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# See also

  - [`PinCFlow.Types.PoissonTypes.Tensor`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Operator`](@ref)

  - [`PinCFlow.Types.PoissonTypes.BiCGSTAB`](@ref)

  - [`PinCFlow.Types.PoissonTypes.Correction`](@ref)
"""
struct Poisson{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: Tensor,
    C <: Operator,
    D <: BiCGSTAB,
    E <: Correction,
}
    rhs::A
    solution::A
    tensor::B
    operator::C
    bicgstab::D
    correction::E
end

function Poisson(domain::Domain)::Poisson
    (; nx, ny, nz) = domain

    (rhs, solution) = (zeros(nx, ny, nz) for i in 1:2)
    tensor = Tensor(domain)
    operator = Operator(domain)
    bicgstab = BiCGSTAB(domain)
    correction = Correction(domain)

    return Poisson(
        rhs,
        solution,
        tensor,
        operator,
        bicgstab,
        correction,
    )
end
