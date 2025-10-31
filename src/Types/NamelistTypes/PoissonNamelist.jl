"""
```julia
PoissonNamelist{A <: Float64, B <: Int, C <: Bool}
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

Construct a `PoissonNamelists` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `tolerance::A`: Tolerance for the convergence criterion of the Poisson solver.

  - `poisson_iterations::B`: Maximum number of iterations performed by the Poisson solver before it terminates regardless of convergence.

  - `preconditioner::C`: Whether to use a preconditioner to accelerate the convergence of the Poisson solver.

  - `dtau::A`: Pseudo-time step coefficient used by the preconditioner.

  - `preconditioner_iterations::B`: Number of iterations performed by the preconditioner.

  - `initial_cleaning::C`: Whether to solve the Poisson problem at initialization to guarantee an initially divergence-free state.

  - `tolerance_is_relative::C`: If set to `true`, the tolerance used for the convergence criterion is given by `tolerance`. If set to `false`, the tolerance is given by `tolerance` divided by a reference value determined from the right-hand side.
"""
struct PoissonNamelist{A <: Float64, B <: Int, C <: Bool}
    tolerance::A
    poisson_iterations::B
    preconditioner::C
    dtau::A
    preconditioner_iterations::B
    initial_cleaning::C
    tolerance_is_relative::C
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
        Float64(tolerance),
        Int(poisson_iterations),
        preconditioner,
        Float64(dtau),
        Int(preconditioner_iterations),
        initial_cleaning,
        tolerance_is_relative,
    )
end
