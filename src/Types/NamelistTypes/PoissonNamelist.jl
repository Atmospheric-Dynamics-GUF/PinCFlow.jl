"""
```julia
PoissonNamelist{
    A <: AbstractFloat,
    B <: Integer,
    C <: Bool,
    D <: AbstractFloat,
    E <: Integer,
    F <: Bool,
    G <: Bool,
}
```

Namelist for parameters used by the Poisson solver.

```julia
PoissonNamelist(;
    tolerance::AbstractFloat = 1.0E-8,
    poisson_iterations::Integer = 1000,
    preconditioner::Bool = true,
    dtau::AbstractFloat = 1.0E+0,
    preconditioner_iterations::Integer = 2,
    initial_cleaning::Bool = true,
    tolerance_is_relative::Bool = false,
)::PoissonNamelist
```

Construct a `PoissonNamelists` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `tolerance::A`: Tolerance for the convergence criterion of the Poisson solver.

  - `poisson_iterations::B`: Maximum number of iterations performed by the Poisson solver before it terminates regardless of convergence.

  - `preconditioner::C`: Whether to use a preconditioner to accelerate the convergence of the Poisson solver.

  - `dtau::D`: Pseudo-time step coefficient used by the preconditioner.

  - `preconditioner_iterations::E`: Number of iterations performed by the preconditioner.

  - `initial_cleaning::F`: Whether to solve the Poisson problem at initialization to guarantee an initially divergence-free state.

  - `tolerance_is_relative::G`: If set to `true`, the tolerance used for the convergence criterion is given by `tolerance`. If set to `false`, the tolerance is given by `tolerance` divided by a reference value determined from the right-hand side.
"""
struct PoissonNamelist{
    A <: AbstractFloat,
    B <: Integer,
    C <: Bool,
    D <: AbstractFloat,
    E <: Integer,
    F <: Bool,
    G <: Bool,
}
    tolerance::A
    poisson_iterations::B
    preconditioner::C
    dtau::D
    preconditioner_iterations::E
    initial_cleaning::F
    tolerance_is_relative::G
end

function PoissonNamelist(;
    tolerance::AbstractFloat = 1.0E-8,
    poisson_iterations::Integer = 1000,
    preconditioner::Bool = true,
    dtau::AbstractFloat = 1.0E+0,
    preconditioner_iterations::Integer = 2,
    initial_cleaning::Bool = true,
    tolerance_is_relative::Bool = false,
)::PoissonNamelist
    return PoissonNamelist(
        tolerance,
        poisson_iterations,
        preconditioner,
        dtau,
        preconditioner_iterations,
        initial_cleaning,
        tolerance_is_relative,
    )
end
