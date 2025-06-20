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
    PoissonNamelist(; <keyword arguments>)

Construct a PoissonNamelist instance, which holds parameters for the Poisson solver.

# Arguments:

  - `tolpoisson::AbstractFloat = 1.0E-8`: Convergence tolerance for the Poisson solver. The solver will terminate when the residual falls below this value.
  - `maxiterpoisson::Integer = 1000`: Maximum number of iterations for the Poisson solver before terminating regardless of convergence.
  - `preconditioner::Bool = true`: Whether to use a preconditioner to accelerate convergence of the Poisson solver.
  - `dtau::AbstractFloat = 1.0`: Time step parameter for the Poisson solver, controls stability and convergence rate.
  - `maxiteradi::Integer = 2`: Maximum number of iterations for the Alternating Direction Implicit (ADI) iterative solver.
  - `initialcleaning::Bool = true`: Whether to perform initial cleaning of the solution field before starting the Poisson solver.
  - `relative_tolerance::Bool = false`: When true, uses relative error for convergence criterion instead of absolute error.
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
