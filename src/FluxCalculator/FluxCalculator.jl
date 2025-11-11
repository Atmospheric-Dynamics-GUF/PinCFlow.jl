"""
```julia
FluxCalculator
```

Module for flux calculation.

Provides functions for MUSCL reconstruction and flux computation.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.Update`](@ref)
"""
module FluxCalculator

using ..Types
using ..Boundaries
using ..Update
using ..PinCFlow

include("apply_1d_muscl!.jl")
include("apply_3d_muscl!.jl")
include("compute_flux.jl")
include("compute_fluxes!.jl")
include("reconstruct!.jl")

export compute_fluxes!, reconstruct!

end
