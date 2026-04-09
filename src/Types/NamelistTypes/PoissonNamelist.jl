"""
```julia
PoissonNamelist
```

Namelist for parameters used by the Poisson solver.

```julia
PoissonNamelist(;
    tolerance::Real = 1.0E-8,
    poisson_iterations::Integer = 1000,
    preconditioner::Bool = true,
    dtau::Real = 1.0E+0,
    preconditioner_iterations::Integer = 2,
    initial_cleaning::Bool = true,
    tolerance_is_relative::Bool = false,
)::PoissonNamelist
```

Construct a `PoissonNamelists` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `tolerance::Float64`: Tolerance for the convergence criterion of the Poisson solver.

  - `poisson_iterations::Int`: Maximum number of iterations performed by the Poisson solver before it terminates regardless of convergence.

  - `preconditioner::Bool`: Whether to use a preconditioner to accelerate the convergence of the Poisson solver.

  - `dtau::Float64`: Pseudo-time step coefficient used by the preconditioner.

  - `preconditioner_iterations::Int`: Number of iterations performed by the preconditioner.

  - `initial_cleaning::Bool`: Whether to solve the Poisson problem at initialization to guarantee an initially divergence-free state.

  - `tolerance_is_relative::Bool`: If set to `true`, the tolerance used for the convergence criterion is given by `tolerance`. If set to `false`, the tolerance is given by `tolerance` divided by a reference value determined from the right-hand side.
"""
struct PoissonNamelist
    tolerance::Float64
    poisson_iterations::Int
    preconditioner::Bool
    dtau::Float64
    preconditioner_iterations::Int
    initial_cleaning::Bool
    tolerance_is_relative::Bool
end

function PoissonNamelist(;
    tolerance::Real = 1.0E-8,
    poisson_iterations::Integer = 1000,
    preconditioner::Bool = true,
    dtau::Real = 1.0E+0,
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
