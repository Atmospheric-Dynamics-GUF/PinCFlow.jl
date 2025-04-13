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
