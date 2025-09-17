"""
```julia
RaySources
```

Module for ray-volume sources.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations`](@ref)
"""
module RaySources

using ...Types
using ..RayOperations
using ..Interpolation # required for activate_multiplewavepackets_source

include("activate_orographic_source!.jl")
include("compute_orographic_mode.jl")
include("activate_multiplewavepackets_source!.jl")

export activate_orographic_source!, activate_multiplewavepackets_source! 

end
