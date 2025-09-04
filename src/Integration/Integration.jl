"""
```julia
Integration
```

Module for integration of the full system.

Provides helper functions for computing the time step, managing time levels and synchronizing fields, as well as the main function for running PinCFlow.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.Update`](@ref)

  - [`PinCFlow.PoissonSolver`](@ref)

  - [`PinCFlow.FluxCalculator`](@ref)

  - [`PinCFlow.Output`](@ref)

  - [`PinCFlow.MSGWaM`](@ref)
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

"""
```julia
Dissipation
```

Singleton for the dissipation step of the turbulence computation.
"""
struct Dissipation end

"""
```julia
Diffusion
```

Singleton for the diffusion step of the turbulence computation.
"""
struct Diffusion end

include("compute_time_step.jl")
include("integrate.jl")
include("modify_compressible_wind!.jl")
include("reset_predictands!.jl")
include("save_backups!.jl")
include("synchronize_compressible_atmosphere!.jl")
include("synchronize_density_fluctuations!.jl")
include("explicit_integration!.jl")
include("implicit_integration!.jl")
include("wkb_integration!.jl")
include("turbulence_computation!.jl")

export Dissipation, Diffusion

export integrate

end
