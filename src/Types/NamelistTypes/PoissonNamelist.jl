"""
```julia
PoissonNamelist{A <: AbstractFloat, B <: Integer, C <: Bool}
```

Namelist for parameters used by the Poisson solver.

```julia
PoissonNamelist(;
    tolpoisson::AbstractFloat = 1.0E-8,
    maxiterpoisson::Integer = 1000,
    preconditioner::Bool = true,
    dtau::AbstractFloat = 1.0E+0,
    maxiteradi::Integer = 2,
    initialcleaning::Bool = true,
    relative_tolerance::Bool = false,
)
```

Construct a `PoissonNamelists` instance with the given keyword arguments as properties.

# Fields/Keywords

- `tolpoisson::A`: Tolerance for the convergence criterion of the Poisson solver.
- `maxiterpoisson::B`: Maximum number of iterations performed by the Poisson solver before it terminates regardless of convergence.
- `preconditioner::C`: Whether to use a preconditioner to accelerate the convergence of the Poisson solver.
- `dtau::A`: Pseudo-time step coefficient used by the preconditioner.
- `maxiteradi::B`: Number of iterations performed by the preconditioner.
- `initialcleaning::C`: Whether to solve the Poisson problem at initialization to guarantee an initially divergence-free state.
- `relative_tolerance::C`: If set to `true`, the tolerance used for the convergence criterion is given by `tolpoisson`. If set to `false`, the tolerance is given by `tolpoisson` divided by a reference value determined from the right-hand side.
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

function PoissonNamelist(;
    tolpoisson::AbstractFloat = 1.0E-8,
    maxiterpoisson::Integer = 1000,
    preconditioner::Bool = true,
    dtau::AbstractFloat = 1.0E+0,
    maxiteradi::Integer = 2,
    initialcleaning::Bool = true,
    relative_tolerance::Bool = false,
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
