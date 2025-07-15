"""
```julia
PoissonNamelist{A <: AbstractFloat, B <: Integer, C <: Bool}
```

Namelist for the Poisson solver (see constructor for parameter descriptions).
"""
struct PoissonNamelist{A <: AbstractFloat, B <: Integer, C <: Bool}
    tolpoisson::A
    maxiterpoisson::B
    preconditioner::C
    dtau::A
    maxiteradi::B
    initialcleaning::C
    relative_tolerance::C
end

"""
```julia
PoissonNamelist(;
    tolpoisson = 1.0E-8,
    maxiterpoisson = 1000,
    preconditioner = true,
    dtau = 1.0E+0,
    maxiteradi = 2,
    initialcleaning = true,
    relative_tolerance = false,
)
```

Construct a PoissonNamelist instance, which holds parameters for the Poisson solver.

# Arguments:

  - `tolpoisson`: Convergence tolerance for the Poisson solver. The solver will terminate when the residual falls below this value.
  - `maxiterpoisson`: Maximum number of iterations for the Poisson solver before terminating regardless of convergence.
  - `preconditioner`: Whether to use a preconditioner to accelerate convergence of the Poisson solver.
  - `dtau`: Time step parameter for the Poisson solver, controls stability and convergence rate.
  - `maxiteradi`: Maximum number of iterations for the Alternating Direction Implicit (ADI) iterative solver.
  - `initialcleaning`: Whether to perform initial cleaning of the solution field before starting the Poisson solver.
  - `relative_tolerance`: When true, uses relative error for convergence criterion instead of absolute error.

# Returns

  - `PoissonNamelist`: Container holding the configuration parameters for the Poisson solver.
"""
function PoissonNamelist(;
    tolpoisson = 1.0E-8,
    maxiterpoisson = 1000,
    preconditioner = true,
    dtau = 1.0E+0,
    maxiteradi = 2,
    initialcleaning = true,
    relative_tolerance = false,
)
    return PoissonNamelist(
        tolpoisson,
        maxiterpoisson,
        preconditioner,
        dtau,
        maxiteradi,
        initialcleaning,
        relative_tolerance,
    )
end
