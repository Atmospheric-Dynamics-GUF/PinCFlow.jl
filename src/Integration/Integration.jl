"""
```julia
Integration
```

Module for integration of the full system.

Provides helper functions for computing the time step, managing time levels and synchronizing fields, as well as the main function for running PinCFlow.
"""
module Integration

using MPI
using Dates
using ..Types
using ..Boundaries
using ..Update
using ..PoissonSolver
using ..FluxCalculator
using ..Output
using ..MSGWaM

include("compute_time_step.jl")
include("integrate.jl")
include("modify_compressible_wind!.jl")
include("reset_predictands!.jl")
include("save_backups!.jl")
include("synchronize_compressible_atmosphere!.jl")
include("synchronize_density_fluctuations!.jl")

export integrate

end
