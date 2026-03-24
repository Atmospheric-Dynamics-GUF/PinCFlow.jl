"""
```julia
BlockedLayer
```

Module for the blocked-layer scheme.

# See also

  - [`PinCFlow.Types`](@ref)

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been validated yet.
"""
module BlockedLayer

using ...Types
using ...PinCFlow

include("compute_blocked_layer!.jl")
include("compute_elevation_difference.jl")
include("compute_slope.jl")
include("include_blocked_flow_drag!.jl")

export compute_blocked_layer!,
    compute_elevation_difference, compute_slope, include_blocked_flow_drag!

end
