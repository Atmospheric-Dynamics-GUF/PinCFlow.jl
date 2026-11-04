"""
```julia
FluxCalculator
```

Module for flux calculation.

Provides functions for MUSCL reconstruction and flux computation.

# Scopes

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Update`](@ref)

# Exports

  - [`PinCFlow.FluxCalculator.compute_fluxes!`](@ref)

  - [`PinCFlow.FluxCalculator.reconstruct!`](@ref)
"""
module FluxCalculator

using ..Macros
using ..Types
using ..Boundaries
using ..Update

include("apply_1d_muscl!.jl")
include("apply_3d_muscl!.jl")
include("compute_flux.jl")
include("compute_fluxes!.jl")
include("reconstruct!.jl")

export compute_fluxes!, reconstruct!

end
