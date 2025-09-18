"""
```julia
PoissonSolver
```

Module for solving the Poisson equation and performing a corrector step.

Provides the functions needed to solve the Poisson equation of the semi-implicit time scheme, using a preconditioned BicGStab algorithm, and correcting the Exner-pressure, momentum and density fluctuations accordingly.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MPIOperations`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.Update`](@ref)
"""
module PoissonSolver

using MPI
using ..Types
using ..MPIOperations
using ..Boundaries
using ..Update
using ..PinCFlow

"""
```julia
Total
```

Singleton for dispatch to application of the full linear operator.
"""
struct Total end

"""
```julia
Horizontal
```

Singleton for dispatch to application of the horizontal part of the linear operator.
"""
struct Horizontal end

include("apply_bicgstab!.jl")
include("apply_corrector!.jl")
include("apply_operator!.jl")
include("apply_preconditioner!.jl")
include("compute_operator!.jl")
include("compute_lhs!.jl")
include("correct!.jl")
include("solve_poisson!.jl")

export apply_corrector!

end
