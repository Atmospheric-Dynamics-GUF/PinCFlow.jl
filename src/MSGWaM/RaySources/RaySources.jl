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

using ..RayOperations
using ...Types
using ...PinCFlow

include("activate_orographic_source!.jl")
include("compute_orographic_mode.jl")
include("activate_continuous_spectral_source!.jl")
include("activate_ray_source!.jl")
export activate_orographic_source!,
      activate_continuous_spectral_source!,
      activate_ray_source!
end
